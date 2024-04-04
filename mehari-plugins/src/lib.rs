/// Contains (placeholders for) structs and traits for the plugin interface.
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct Dummy {
    pub some_field: String,
}

impl Dummy {
    pub fn change_something(&mut self, dummy: &str) {
        self.some_field = dummy.into();
    }
}
