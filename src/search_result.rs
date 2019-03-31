use std::cmp::Ordering;

#[derive(PartialEq, Debug, Clone)]
pub struct SearchResult {
    pub node: usize,
    pub distance: f32
}

impl SearchResult {

    pub fn new(node: usize, distance: f32) -> SearchResult {
        SearchResult{ node, distance }
    }

    pub fn none() -> SearchResult {
        SearchResult::new(0, std::f32::INFINITY)
    }

}

impl Eq for SearchResult {}

impl PartialOrd for SearchResult {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.distance < other.distance {
            Some(Ordering::Less)
        } else if self.distance > other.distance {
            Some(Ordering::Greater)
        }
        else {
            None
        }
    }
}

impl Ord for SearchResult {
    fn cmp(&self, other: &SearchResult) -> Ordering {
        self.partial_cmp(other).unwrap_or_else(|| Ordering::Equal)
    }
}