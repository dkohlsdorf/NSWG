use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;

pub struct UCRInstance {
    pub label: usize,
    pub instance: Vec<f32>
}

impl UCRInstance {

    pub fn from_file(file: String) -> (Vec<UCRInstance>, usize, usize) {
        let mut max_label = 0;
        let f = File::open(file).unwrap();
        let reader = BufReader::new(f);
        let mut dataset = vec![];
        let mut dim = 0;
        for line in reader.lines() {
            let line = line.unwrap();
            let mut label = 0;
            let mut instance = vec![];
            for (i, cmp) in line.split(",").enumerate() {
                if i == 0 {
                    label = cmp.parse::<usize>().unwrap();
                    if label > max_label {
                        max_label = label;
                    }
                } else {
                    instance.push(cmp.parse::<f32>().unwrap());
                }
            }
            dim = instance.len();
            dataset.push(UCRInstance {label: label, instance: instance});
        }
        (dataset, max_label + 1, dim)
    }

}