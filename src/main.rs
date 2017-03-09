extern crate rand;
use rand::Rng;
use rand::distributions::{IndependentSample, Range};

#[derive(Clone)]
struct Population<T> {
    individuals: Vec<Individual<T>>,
    size: i32,
    n_gens: i32,
    ub: T,
    lb: T,
}

#[derive(Clone)]
struct Individual<T> {
    size: i32,
    gene: Vec<T>,
    ub: T,
    lb: T,
}

impl<T> Population<T> {
    fn new (size: i32, n_gens:i32, lb:T, ub:T) -> Population<T> {
        Population::<T> {
            individuals: Vec::new(),
            size: size,
            n_gens: n_gens,
            lb: lb,
            ub: ub,
        }
    }
}

impl Population<i32> {
    fn init(&mut self) {
        for _ in 0..self.size {
            let mut w:Individual<i32> = Individual::<i32>::new(self.n_gens, self.lb, self.ub);
            w.birth();
            //w.print();
            self.individuals.push(w);
        }
    }
}

impl Population<bool> {
    fn init(&mut self) {
        for _ in 0..self.size {
            let mut w:Individual<bool> = Individual::<bool>::new(self.n_gens, self.lb, self.ub);
            w.birth();
            //w.print();
            self.individuals.push(w);
        }
    }
}

impl Population<f64> {
    fn init(&mut self) {
        for _ in 0..self.size {
            let mut w:Individual<f64> = Individual::<f64>::new(self.n_gens, self.lb, self.ub);
            w.birth();
            //w.print();
            self.individuals.push(w);
        }
    }
}

impl Individual<i32> {
    fn birth(&mut self) {
        let between = Range::new(self.lb, self.ub);
        let mut rng = rand::thread_rng();

        for _ in 0..self.size {
            let w: i32 = between.ind_sample(&mut rng);
            self.gene.push(w);
        }
    }

    fn print(&self) {
        for i in 0..self.gene.len() {
            println!("{0:4} {1}", i, self.gene[i]);
        }
    }
}

impl Individual<f64> {
    fn birth(&mut self) {
        let between = Range::new(self.lb, self.ub);
        let mut rng = rand::thread_rng();

        for _ in 0..self.size {
            let w: f64 = between.ind_sample(&mut rng);
            self.gene.push(w);
        }
    }

    fn print(&self) {
        for i in 0..self.gene.len() {
            println!("{0:4} {1}", i, self.gene[i]);
        }
    }
}

impl Individual<bool> {
    fn birth(&mut self) {
        let mut rng = rand::thread_rng();

        for _ in 0..self.size {
            let w: bool = rng.gen::<bool>();
            self.gene.push(w);
        }
    }

    fn print(&self) {
        for i in 0..self.gene.len() {
            println!("{0:4} {1}", i, self.gene[i]);
        }
    }
}

impl<T> Individual<T> {
    fn new(size:i32, lb:T, ub:T) -> Individual<T> {
        Individual::<T> {
            size: size,
            gene: Vec::new(),
            ub: ub,
            lb: lb,
        }
    }
}

fn main () {
    let mut pop:Population<f64> = Population::<f64>::new(10, 10, -5.12, 5.12);
    pop.init();

    for x in &pop.individuals {
        x.print();
    }
}
