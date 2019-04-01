pub mod nswg;
pub mod result_store;
pub mod search_result;
pub mod ucr_ts;
use std::time::Instant;

pub fn knn(neighbors: Vec<usize>, labels: usize) -> usize {
    let mut labels = vec![0; labels];
    for neighbor in neighbors {
        labels[neighbor] += 1;
    }
    let mut max_label = 0;
    let mut max_value = 0;
    for (l, c) in labels.iter().enumerate() {
        if *c > max_value {
            max_value = *c;
            max_label = l;
        }
    }
    max_label
}

fn main() {
    println!("Navigable Small World Graphs");
    let (train, max_labels, dim) = ucr_ts::UCRInstance::from_file(String::from(
        "PATH TO TRAIN"
    ));
    let (test, _, _) = ucr_ts::UCRInstance::from_file(String::from(
        "PATH TO TEST",
    ));
    let mut nswg = nswg::NavigableSmallWorldGraph::new(dim);
    let now = Instant::now();
    for i in 0..train.len() {
        if i % 100 == 0 && i > 0 {
            println!(
                "\tINSERT: {} / {} {} [ms] {} [ms]",
                i,
                train.len(),
                now.elapsed().as_millis(),
                now.elapsed().as_millis() / i as u128
            );
        }
        nswg.insert(i, &train[i].instance[..], 25, 35, None);
    }
    println!(
        "Done instering {} {} [ms] {} dim ",
        train.len(),
        now.elapsed().as_millis(),
        dim
    );
    let mut n_correct = 0.0;
    let mut total_steps = 0;
    for (i, x) in test.iter().enumerate() {
        let (neighbors, n_steps) = nswg.search(&x.instance[..], 25, 1, None);
        total_steps += n_steps;
        if i % 100 == 0 && i > 0 {
            println!(
                "\tSEARCH: {} / {} ... steps / naive {} {}",
                i,
                test.len(),
                total_steps / i,
                train.len()
            );
        }
        let labels: Vec<usize> = neighbors.iter().map(|n| train[n.node].label).collect();
        if x.label == knn(labels, max_labels) {
            n_correct += 1.0;
        }
    }
    n_correct /= test.len() as f32;
    println!(
        "{} [ms] with {} error",
        now.elapsed().as_millis(),
        1.0 - n_correct
    );
}
