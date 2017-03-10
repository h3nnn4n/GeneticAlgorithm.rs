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
    fitness: f64,
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
        for ref mut x in &mut self.individuals {
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

    fn mutate_everyone(&mut self){
        for ref mut x in &mut self.individuals {
            x.mutate();
        }
    }

    fn print_pop_fit(&self){
        for ref x in &self.individuals {
            x.print_fitness();
        }
        println!();
    }

    fn print_best_fit(&self){
        let mut best = self.individuals[0].fitness;
        let mut best_obj = self.individuals[0].objective_function();
        for ref x in &self.individuals {
            let w = x.fitness;
            if w > best {
                best_obj = x.objective_function();
                best = w;
            }
        }
        println!("{} {}", best, best_obj );
        //println!();
    }

    fn get_best(&self) -> Individual<f64> {
        let mut pos = 0;
        for i in 1..self.individuals.len(){
            if self.individuals[pos].fitness < self.individuals[i].fitness {
                pos = i;
            }
        }

        return self.individuals[pos].clone();
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

        self.fitness = self.fitness();
    }

    fn objective_function(&self) -> f64 {
        return rastrigin(self.gene.as_slice());
    }

    fn fitness(&self) -> f64 {
        let f = self.objective_function();
        if f < 0.0001 {
            return 100000.;
        } else {
            return 1.0/f;
        }
    }

    fn print(&self) {
        for i in 0..self.gene.len() {
            println!("{0:4} {1}", i, self.gene[i]);
        }
        self.print_fitness();
    }

    fn print_fitness(&self) {
        println!("Fitness: {}", self.fitness());
    }

    fn mutate(&mut self){
        //let p = rand::thread_rng().gen_range(0, self.gene.len());
        let mut changed:bool = false;
        for p in 0..self.gene.len() {
            if rand::thread_rng().gen::<f64>() < self.mutation_chance {
                self.gene[p] = self.gene[p] + rand::thread_rng().gen::<f64>() - 0.5;
                changed = true;
            }
        }

        if changed {
            self.fitness = self.fitness();
        }
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
            fitness: -1.,
        }
    }

    fn set_mutation_chance(&mut self, chance: f64){
        self.mutation_chance = chance;
    }
}

fn main () {
    let n_dim = 2;
    let pop_size = 10;
    let max_iter = 1000;

    let mut pop:Population<f64> = Population::<f64>::new(pop_size, n_dim, -5.12, 5.12);

    pop.init();
    pop.set_mutation_chance(0.1);

    //pop.print_pop_fit();

    for _ in 0..max_iter{
        let best = pop.get_best();
        pop.print_best_fit();
        pop.mutate_everyone();

        pop.individuals[0] = best;
    }
}
