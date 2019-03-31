use std::collections::BTreeSet;
use crate::search_result::*;

pub struct ResultStore {
    ordered: BTreeSet<SearchResult>
}

impl ResultStore {    

    pub fn new() -> ResultStore {
        ResultStore { ordered: BTreeSet::new() }
    }

    pub fn len(&self) -> usize {
        self.ordered.len()
    }

    pub fn take(&self, k: usize) -> Vec<SearchResult> {
        self.ordered.iter().take(k).map(|x| x.clone()).collect()
    }

    pub fn best(&self, at: usize) -> SearchResult {
        let k_best = self.ordered.iter().take(at);
        if k_best.len() < at {
            SearchResult::none()
        } else {
            k_best.last().unwrap().clone()
        }
    }

    pub fn insert(&mut self, result: SearchResult) {
        self.ordered.insert(result);
    }
    
    pub fn remove(&mut self, result: &SearchResult) {
        self.ordered.remove(result);
    }

}