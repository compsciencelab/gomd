package sim

import (
	"fmt"
	"io/ioutil"
	"math"
	"os"
	"strconv"
	"strings"
)

const (
	TimeFactor Real = 48.88821
	Boltzman   Real = 0.001987191
)

type Real float32

type Vec3 struct {
	X, Y, Z Real
}

type Sim struct {
	pos, vel, force    []Vec3
	natoms             int
	dt, mass, Kin, Pot Real
	sigma, eps         Real
	box                Vec3
}

func New(natoms int, box Vec3) (sim *Sim) {
	//var sim Sim
	sim = new(Sim)
	sim.pos = make([]Vec3, natoms, natoms)
	sim.vel = make([]Vec3, natoms, natoms)
	sim.force = make([]Vec3, natoms, natoms)
	sim.Kin = 0
	sim.Pot = 0
	sim.natoms = natoms
	sim.sigma = Real(3.4) //do I need an explicit conversion as the type is already defined?
	sim.eps = Real(0.238)
	sim.dt = 1.0 / TimeFactor
	sim.mass = Real(39.948)
	sim.box = box
	return sim
}

func (sim *Sim) Read(filename string) { //this is to read a tabular file. I start to love Python
	content, err := ioutil.ReadFile(filename)
	if err != nil {
		panic(err)
	}
	result := strings.Split(string(content), "\n") //splits seems to add and extra line, so I check num values
	for i := range result {
		values := strings.Fields(result[i])
		if len(values) == 3 {
			v0, _ := strconv.ParseFloat(values[0], 64)
			v1, _ := strconv.ParseFloat(values[1], 64)
			v2, _ := strconv.ParseFloat(values[2], 64)
			sim.pos[i].X = Real(v0)
			sim.pos[i].Y = Real(v1)
			sim.pos[i].Z = Real(v2)
		}
	}

}

func (sim *Sim) Write(filename string) {
	fo, err := os.Create(filename)
	if err != nil {
		panic(err)
	}
	for i := range sim.pos {
		fmt.Fprintf(fo, "%f %f %f\n", sim.pos[i].X, sim.pos[i].Y, sim.pos[i].Z)
	}

}

func (sim *Sim) GetTemperature() Real {
	return sim.Kin * 2.0 / (3.0 * Boltzman * Real(sim.natoms))
}

func (sim *Sim) FirstVV() {
	//First part of Velocity Verlet
	imass := 1.0 / sim.mass
	dt := sim.dt
	//pos := &sim.pos  //can I do name aliasing?
	for i := 0; i < sim.natoms; i++ {
		coeff := 0.5 * dt * dt * imass
		sim.pos[i].X += dt*sim.vel[i].X + coeff*sim.force[i].X
		sim.pos[i].Y += dt*sim.vel[i].Y + coeff*sim.force[i].Y
		sim.pos[i].Z += dt*sim.vel[i].Z + coeff*sim.force[i].Z

		coeffvel := 0.5 * dt * imass
		sim.vel[i].X += coeffvel * sim.force[i].X
		sim.vel[i].Y += coeffvel * sim.force[i].Y
		sim.vel[i].Z += coeffvel * sim.force[i].Z
	}
}

func (sim *Sim) SecondVV() {
	imass := Real(1.0 / sim.mass) //is this Real or double?
	tmpE := Real(0.0)
	coeffvel := Real(0.5) * sim.dt * imass
	for i := 0; i < sim.natoms; i++ {
		//vec.Inc(sim.vel[i],vec.Scale(coeffvel,sim.force[i]))
		sim.vel[i].X += coeffvel * sim.force[i].X
		sim.vel[i].Y += coeffvel * sim.force[i].Y
		sim.vel[i].Z += coeffvel * sim.force[i].Z
		tmpE += sim.vel[i].X*sim.vel[i].X + sim.vel[i].Y*sim.vel[i].Y + sim.vel[i].Z*sim.vel[i].Z
	}
	sim.Kin = tmpE * 0.5 * sim.mass
}

func (sim *Sim) ResetForce() {
	for i := 0; i < sim.natoms; i++ {
		sim.force[i] = Vec3{0, 0, 0}
	}

}

func round(x Real) Real {
	return Real(math.Floor(float64(x)))
}

func dist2pbc(dist, box Vec3) Real {
	v := Vec3{
		dist.X - box.X*round(dist.X/box.X),
		dist.Y - box.Y*round(dist.Y/box.Y),
		dist.Y - box.Y*round(dist.Y/box.Y)}
	return v.X*v.X + v.Y*v.Y + v.Z*v.Z
}

func (sim *Sim) ComputeNonBonded() {
	Epot := Real(0)
	rcut2 := Real(12.0 * 12.0)
	A := sim.eps * Real(4.0*math.Pow(float64(sim.sigma), 12.0))
	B := sim.eps * Real(4.0*math.Pow(float64(sim.sigma), 6.0))
	for i := 0; i < sim.natoms; i++ {
		for j := i + 1; j < sim.natoms; j++ {
			rij := Vec3{sim.pos[i].X - sim.pos[j].X, sim.pos[i].Y - sim.pos[j].Y, sim.pos[i].Z - sim.pos[j].Z}
			//r2 := rij.X*rij.X + rij.Y*rij.Y + rij.Z*rij.Z
			r2 := dist2pbc(rij, sim.box)
			if r2 < rcut2 {
				rminus1 := 1.0 / Real(math.Sqrt(float64(r2)))
				rminus2 := rminus1 * rminus1
				rminus6 := rminus2 * rminus2 * rminus2
				rminus12 := rminus6 * rminus6
				AmbTerm := (A*rminus6 - B) * rminus6
				Epot += AmbTerm
				forcer := (Real(6.0)*(A*rminus12+AmbTerm)*rminus1 - AmbTerm) * rminus1
				forceij := Vec3{forcer * rij.X, forcer * rij.Y, forcer * rij.Z}
				sim.force[i].X += forceij.X
				sim.force[j].X -= forceij.X
				sim.force[i].Y += forceij.Y
				sim.force[j].Y -= forceij.Y
				sim.force[i].Z += forceij.Z
				sim.force[j].Z -= forceij.Z
			}
		}
	}
	sim.Pot = Epot
}
