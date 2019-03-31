extern crate rand;

use rand::Rng;
use std::collections::{HashMap, HashSet};

use crate::search_result::*;
use crate::result_store::*;

/// A navigable small world graph defined
/// by a sparse edge set and a flat vector set, one for each vertex.
pub struct NavigableSmallWorldGraph {
    pub edges: HashMap<usize, Vec<usize>>,
    pub instances: Vec<f32>,
    pub dim: usize
}

impl NavigableSmallWorldGraph {

    /// Create an empty graph
    /// * `dim` dimension of the node's vectors
    pub fn new(dim: usize) -> NavigableSmallWorldGraph {
        NavigableSmallWorldGraph {
            edges: HashMap::new(),
            instances: Vec::new(),
            dim: dim        
        }
    }

    /// number of instances in the graph
    pub fn n_instances(&self) -> usize {
        self.instances.len() / self.dim
    }

    /// Insert a vector into the graph, the nearest neihbors to the vector become the edges.
    /// * `id` vector id
    /// * `query` the vector
    /// * `n_searches` number of search restarts
    /// * `k_neighbors` number of search results, here number of edges
    /// * `align_band` Sakoe Chiba Band for Dynamic Time Warping, if not set we use Euclidean Distance
    pub fn insert(&mut self, id: usize, query: &[f32], n_searches: usize, k_neighbors: usize, align_band: Option<usize>) {
        if self.n_instances() > k_neighbors {
            let knn = self.search(query, n_searches, k_neighbors, align_band).0
                .iter()
                .map(|result| result.node)
                .collect();
            self.set(id, knn, query);
        } else {
            let knn: Vec<usize> = (0 .. self.n_instances()).collect();
            self.set(id, knn, query);
        }
    }

    /// Given a query find the k-nearest neighbors using the nswg graph
    /// * `query` The query vector
    /// * `n_searches` number of restarts
    /// * `k_neighbors` number of results
    /// * `align_band` Sakoe Chiba Band for Dynamic Time Warping, if not set we use Euclidean Distance
    pub fn search(&self, query: &[f32], n_searches: usize, k_neighbors: usize, align_band: Option<usize>) -> (Vec<SearchResult>, usize){
        let mut candidates = ResultStore::new();
        let mut results    = ResultStore::new();
        let mut visited    = HashSet::new();
        let mut rng        = rand::thread_rng();   
        let mut n_steps    = 0;     
        for _attempt in 0 .. n_searches {
            let mut bsf     = ResultStore::new();
            
            // stop if we visited all of the nodes
            if self.n_instances() == visited.len() {
                break;
            }

            // sample entry point, make sure we never used it before
            // and insert the point as a candidate
            let mut entry_point = rng.gen_range(0, self.n_instances());
            while visited.contains(&entry_point) {
                entry_point = rng.gen_range(0, self.n_instances());
            }
            let distance    = self.distance(query, entry_point, align_band);         
            candidates.insert(SearchResult::new(entry_point, distance));

            while candidates.len() > 0 {
                // check if we can stop: abendon search if the result is the worst than the top k so far
                let candidate = candidates.best(1);
                let kth = f32::min(results.best(k_neighbors).distance, bsf.best(k_neighbors).distance);
                candidates.remove(&candidate);
                if candidate.distance > kth {
                    break;
                }
                // add expand all edges 
                for edge in self.edges[&candidate.node].iter() {                                        
                    if !visited.contains(edge) {
                        n_steps += 1;
                        visited.insert(edge);
                        let edge     = *edge;
                        let distance = self.distance(query, edge, align_band);                                                                        
                        candidates.insert(SearchResult::new(edge, distance));
                        bsf.insert(SearchResult::new(edge, distance));                    
                    }
                }                
            }           
            for res in bsf.take(k_neighbors) {
                results.insert(res);
            }
        }
        (results.take(k_neighbors), n_steps)
    }

    fn set(&mut self, id: usize, neighbors: Vec<usize>, instance: &[f32]) {
        self.edges.insert(id, Vec::new());
        for neighbor in neighbors {
            self.edges.get_mut(&id).unwrap().push(neighbor);
            self.edges.get_mut(&neighbor).unwrap().push(id);
        }
        self.instances.extend_from_slice(instance);
    }

    fn distance(&self, query: &[f32], node: usize, align_band: Option<usize>) -> f32 {
        match align_band {
            Some(w) => self.dtw(query, node, w),
            None    => self.euclidean(query, node)
        }
    }

    fn euclidean(&self, query: &[f32], node: usize) -> f32 {
        assert!(query.len() == self.dim);
        let y = &self.instances[node  * self.dim .. (node + 1) * self.dim];
        let mut distance = 0.0;
        for i in 0 .. query.len() {
            distance += f32::powf(query[i] - y[i], 2.0);
        }
        f32::sqrt(distance)
    }

    fn dtw(&self, query: &[f32], node: usize, w: usize) -> f32 {
        let y  = &self.instances[node  * self.dim .. (node + 1) * self.dim];
        let n  = query.len();
        let m  = y.len();
        // allocating here is the main bottleneck
        let mut dp = vec![std::f32::INFINITY; (n + 1) * (m + 1)]; 
        dp[0] = 0.0;
        for i in 1 .. n + 1 {
            for j in usize::max(NavigableSmallWorldGraph::sub(i, w), 1) .. usize::min(i + w, m + 1) {
                let distance = f32::powf(query[i - 1] - y[j - 1], 2.0);                
                dp[i * (n + 1) + j] = distance + NavigableSmallWorldGraph::min3(
                    dp[(i - 1) * (n + 1) + j],
                    dp[i       * (n + 1) + (j - 1)], 
                    dp[(i - 1) * (n + 1) + (j - 1)]
                )
            }
        }
        dp[dp.len() - 1]
    }
    
    fn sub(x: usize, y: usize) -> usize{
        if x > y {
            x - y
        } else {
            0
        }
    }

    fn min3(x: f32, y: f32, z: f32) -> f32 {
        let mut min = x;
        if y < min {
            min = y;
        }
        if z < min {
            min = z;
        }
        min
    }

}
