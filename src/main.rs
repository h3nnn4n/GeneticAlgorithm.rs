extern crate rand;
use rand::Rng;
use rand::distributions::{IndependentSample, Range};
use std::f64::consts::PI;

#[allow(dead_code)]
fn rastrigin(pos: &[f64]) -> f64 {
    return pos.iter().fold(0.0, |sum, x| sum + x * x - 10.0 * (2.0 * PI * x).cos() + 10.0);
}

#[derive(Clone)]
struct Population<T> {
    individuals: Vec<Individual<T>>,
    size: i32,
    n_gens: i32,
    ub: T,
    lb: T,
    mutation_chance: f64,
    crossover_chance: f64,
}

#[derive(Clone)]
struct Individual<T> {
    size: i32,
    gene: Vec<T>,
    ub: T,
    lb: T,
    mutation_chance: f64,
}

#[allow(dead_code)]
impl<T> Population<T> {
    fn new (size: i32, n_gens:i32, lb:T, ub:T) -> Population<T> {
        Population::<T> {
            individuals: Vec::new(),
            size: size,
            n_gens: n_gens,
            lb: lb,
            ub: ub,
            mutation_chance: 0.05,
            crossover_chance: 0.8,
        }
    }

    fn set_mutation_chance(&mut self, chance: f64){
        self.mutation_chance = chance;
        for &mut x in self.individuals {
            x.set_mutation_chance(chance);
        }
    }
}

#[allow(dead_code)]
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

#[allow(dead_code)]
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

#[allow(dead_code)]
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

#[allow(dead_code)]
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

#[allow(dead_code)]
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

    fn mutate(&mut self){
        unimplemented!();
    }

    fn fitness(self) -> f64 {
        return -rastrigin(self.gene.as_slice());
    }
}

#[allow(dead_code)]
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

#[allow(dead_code)]
impl<T> Individual<T> {
    fn new(size:i32, lb:T, ub:T) -> Individual<T> {
        Individual::<T> {
            size: size,
            gene: Vec::new(),
            ub: ub,
            lb: lb,
            mutation_chance: 0.05,
        }
    }

    fn set_mutation_chance(&mut self, chance: f64){
        self.mutation_chance = chance;
    }
}

fn main () {
    let n_dim = 10;
    let pop_size = 10;
    let mut pop:Population<f64> = Population::<f64>::new(pop_size, n_dim, -5.12, 5.12);
    pop.init();
    pop.set_mutation_chance(0.1);

    for mut x in pop.individuals {
        x.print();
        println!();
    }
}
