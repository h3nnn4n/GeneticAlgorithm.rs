extern crate rand;
use rand::distributions::{IndependentSample, Range};

#[derive(Clone)]
struct Individual<T> {
    size: i32,
    gene: Vec<T>,
    ub: T,
    lb: T,
}

impl Individual<f64> {
    fn init (&mut self, size:i32) {
        self.size = size;
    }

    fn birth(&mut self) {
        let between = Range::new(self.lb, self.ub);
        let mut rng = rand::thread_rng();

        for i in 0..self.size {
            let w: f64 = between.ind_sample(&mut rng);
        }
    }
}

fn main () {
    let a: i32 = 1;
    let b = a;
}
