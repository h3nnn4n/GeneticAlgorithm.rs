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

    fn init(&mut self) {
        for _ in 0..self.size {
            let w:Individual<T> = Individual::<T>::new(&self.n_gens, &mut self.lb, &mut self.ub);
            //w.birth();
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
    fn new(size:&i32, lb:&mut T, ub:&mut T) -> Individual<T> {
        Individual::<T> {
            size: *size,
            gene: Vec::new(),
            ub: *ub,
            lb: *lb,
        }
    }
}

fn main () {
    //let mut ind:Individual<bool> = Individual::new(10, false, true);
    //ind.birth();
    //ind.print();
}
