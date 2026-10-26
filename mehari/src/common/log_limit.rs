//! Log the first warnings and errors of each call site, and count the rest.

use std::fmt;
use std::sync::{Arc, Mutex, MutexGuard, PoisonError};

use indexmap::IndexMap;
use thousands::Separable;
use tracing::callsite::Identifier;
use tracing::field::{Field, Visit};
use tracing::{Event, Level, Subscriber};
use tracing_subscriber::layer::{Context, Layer};

/// Warnings and errors per call site that are logged one by one.
const LOGGED_PER_CALLSITE: usize = 5;

/// The target of the summary, which is never limited.
const SUMMARY_TARGET: &str = "mehari::log_limit";

/// A layer that drops the warnings and errors of a call site after the first
/// `LOGGED_PER_CALLSITE`, for all layers; `log_summary` reports what it dropped.
#[derive(Clone, Default)]
pub struct LogLimit {
    seen: Arc<Mutex<IndexMap<Identifier, Seen>>>,
}

struct Seen {
    count: usize,
    level: Level,
    /// The message of the first event.
    first: String,
}

impl<S: Subscriber> Layer<S> for LogLimit {
    fn event_enabled(&self, event: &Event<'_>, _ctx: Context<'_, S>) -> bool {
        let meta = event.metadata();
        if *meta.level() > Level::WARN || meta.target() == SUMMARY_TARGET {
            return true;
        }
        let count = (self.seen().get_mut(&meta.callsite())).map(|seen| {
            seen.count += 1;
            seen.count
        });
        let count = count.unwrap_or_else(|| {
            // The fields can log while they are formatted, which takes the lock again.
            let first = message(event);
            let mut seen = self.seen();
            let seen = seen.entry(meta.callsite()).or_insert_with(|| Seen {
                count: 0,
                level: *meta.level(),
                first,
            });
            seen.count += 1;
            seen.count
        });
        count <= LOGGED_PER_CALLSITE
    }
}

impl LogLimit {
    fn seen(&self) -> MutexGuard<'_, IndexMap<Identifier, Seen>> {
        self.seen.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// Log, per call site and at its level, how many events were dropped, and forget them.
    pub fn log_summary(&self) {
        let seen = std::mem::take(&mut *self.seen());
        for seen in seen.into_values() {
            if seen.count <= LOGGED_PER_CALLSITE {
                continue;
            }
            let more = (seen.count - LOGGED_PER_CALLSITE).separate_with_commas();
            if seen.level == Level::ERROR {
                tracing::error!(target: SUMMARY_TARGET, "{more} more like: {}", seen.first);
            } else {
                tracing::warn!(target: SUMMARY_TARGET, "{more} more like: {}", seen.first);
            }
        }
    }
}

/// Await `command`, then log the summary of `log_limit`, also if `command` failed.
pub async fn summarize<T>(log_limit: Option<&LogLimit>, command: impl Future<Output = T>) -> T {
    let result = command.await;
    if let Some(log_limit) = log_limit {
        log_limit.log_summary();
    }
    result
}

/// The `message` field of `event`.
fn message(event: &Event<'_>) -> String {
    struct Message(String);
    impl Visit for Message {
        fn record_debug(&mut self, field: &Field, value: &dyn fmt::Debug) {
            if field.name() == "message" {
                self.0 = format!("{value:?}");
            }
        }
    }
    let mut visitor = Message(String::new());
    event.record(&mut visitor);
    visitor.0
}

#[cfg(test)]
mod test {
    use std::sync::{Arc, Mutex};

    use tracing::{Event, Level, Subscriber};
    use tracing_subscriber::layer::{Context, Layer, SubscriberExt};

    use super::LogLimit;

    /// The level and message of the events that reach it.
    #[derive(Clone, Default)]
    struct Logged(Arc<Mutex<Vec<(Level, String)>>>);

    impl<S: Subscriber> Layer<S> for Logged {
        fn on_event(&self, event: &Event<'_>, _ctx: Context<'_, S>) {
            let mut logged = self.0.lock().unwrap();
            logged.push((*event.metadata().level(), super::message(event)));
        }
    }

    /// Warnings and errors are dropped per call site after the first 5, other levels are not.
    /// The summary counts the dropped events at their level.
    #[test]
    fn limit_per_call_site() {
        let log_limit = LogLimit::default();
        let logged = Logged::default();
        let subscriber = tracing_subscriber::registry()
            .with(log_limit.clone())
            .with(logged.clone());

        tracing::subscriber::with_default(subscriber, || {
            for i in 0..6 {
                tracing::warn!("warning {i}");
            }
            for i in 0..6 {
                tracing::info!("info {i}");
            }
            for i in 0..1007 {
                tracing::error!("error {i}");
            }
            log_limit.log_summary();
        });

        let expected: Vec<_> = (0..5)
            .map(|i| (Level::WARN, format!("warning {i}")))
            .chain((0..6).map(|i| (Level::INFO, format!("info {i}"))))
            .chain((0..5).map(|i| (Level::ERROR, format!("error {i}"))))
            .chain([
                (Level::WARN, "1 more like: warning 0".into()),
                (Level::ERROR, "1,002 more like: error 0".into()),
            ])
            .collect();
        assert_eq!(*logged.0.lock().unwrap(), expected);
    }

    /// The summary comes also after a command that fails.
    #[test]
    fn summary_after_failure() {
        let log_limit = LogLimit::default();
        let logged = Logged::default();
        let subscriber = tracing_subscriber::registry()
            .with(log_limit.clone())
            .with(logged.clone());

        let result = tracing::subscriber::with_default(subscriber, || {
            futures::executor::block_on(super::summarize(Some(&log_limit), async {
                for i in 0..6 {
                    tracing::warn!("warning {i}");
                }
                Err::<(), _>("failed")
            }))
        });

        assert_eq!(result, Err("failed"));
        assert_eq!(
            logged.0.lock().unwrap().last(),
            Some(&(Level::WARN, "1 more like: warning 0".to_string()))
        );
    }
}
