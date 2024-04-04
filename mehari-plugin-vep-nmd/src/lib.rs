use extism_pdk::*;
use mehari_plugins::*;

#[plugin_fn]
pub fn process(record: String) -> FnResult<String> {
    let mut record: Dummy = serde_json::from_str(&record).expect("Failed deserializing record");
    record.change_something("after");
    let record = serde_json::to_string(&record).expect("Failed serializing record");
    Ok(record)
}
