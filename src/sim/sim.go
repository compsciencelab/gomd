package sim

import "math"

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
}

func New(natoms int) (sim *Sim) {
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
	return sim
}

func Read(sim *Sim) {

}

func Write(sim *Sim) {

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
func (sim *Sim) ComputeNonBonded() {
	Epot := Real(0)
	rcut2 := Real(12.0 * 12.0)
	A := sim.eps * Real(4.0*math.Pow(float64(sim.sigma), 12.0))
	B := sim.eps * Real(4.0*math.Pow(float64(sim.sigma), 6.0))
	for i := 0; i < sim.natoms; i++ {
		for j := i + 1; j < sim.natoms; j++ {
			var rij Vec3
			rij.X = sim.pos[i].X - sim.pos[j].X
			rij.Y = sim.pos[i].Y - sim.pos[j].Y
			rij.Z = sim.pos[i].Z - sim.pos[j].Z
			r2 := rij.X*rij.X + rij.Y*rij.Y + rij.Z*rij.Z
			if r2 < rcut2 {
				r := Real(math.Sqrt(float64(r2)))
				r_1 := 1.0 / r
				r_2 := r_1 * r_1
				r_6 := r_2 * r_2 * r_2
				r_12 := r_6 * r_6

				AmbTerm := (A*r_6 - B) * r_6
				Epot += AmbTerm
				force_r := (Real(6.0)*(A*r_12+AmbTerm)*r_1 - AmbTerm) * r_1

				var forceij Vec3
				forceij.X = force_r * rij.X
				sim.force[i].X += forceij.X
				sim.force[j].X -= forceij.X
				forceij.Y = force_r * rij.Y
				sim.force[i].Y += forceij.Y
				sim.force[j].Y -= forceij.Y
				forceij.Z = force_r * rij.Z
				sim.force[i].Z += forceij.Z
				sim.force[j].Z -= forceij.Z
			}
		}
	}
	sim.Pot = Epot
}
