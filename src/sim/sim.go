package sim

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

func New(natoms int) Sim {
	var sim Sim
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
	coeffvel := 0.5 * sim.dt * imass
	//Velocity update, second part//
	for i := 0; i < sim.natoms; i++ {
		sim.vel[i].X += coeffvel * sim.force[i].X
		sim.vel[i].Y += coeffvel * sim.force[i].Y
		sim.vel[i].Z += coeffvel * sim.force[i].Z
		tmpE += sim.vel[i].X*sim.vel[i].X + sim.vel[i].Y*sim.vel[i].Y + sim.vel[i].Z*sim.vel[i].Z
	}
	sim.Kin = tmpE * 0.5 * sim.mass
}

/*
func ComputeNonBonded(sim *Sim) {
	Epot := Real(0)
	rcut2 := Real(12.0 * 12.0)
	pos := &sim.pos     //correcto?
	force := &sim.force //correcto?
	A := Real(4.0 * sim.eps * Pow(sim.sigma, 12.0))
	B := Real(4.0 * sim.eps * Pow(sim.sigma, 6.0))
	for i := 0; i < sim.natoms; i++ {
		for j := i + 1; j < natoms; j++ {
			var rij Vec3
			rij.x = pos[i].x - pos[j].x
			rij.y = pos[i].y - pos[j].y
			rij.z = pos[i].z - pos[j].z
			r2 := rij.x*rij.x + rij.y*rij.y + rij.z*rij.z
			if r2 < rcut2 {
				r := Real(sqrt(r2))
				r_1 := 1.0 / r
				r_2 := r_1 * r_1
				r_6 := r_2 * r_2 * r_2
				r_12 := r_6 * r_6

				AmbTerm := (A*r_6 - B) * r_6
				Epot += AmbTerm
				force_r = (6.0*(A*r_12+AmbTerm)*r_1 - AmbTerm) * r_1

				var forceij Vec3
				forceij.x = force_r * rij.x
				force[i].x += forceij.x
				force[j].x -= forceij.x
				forceij.y = force_r * rij.y
				force[i].y += forceij.y
				force[j].y -= forceij.y
				forceij.z = force_r * rij.z
				force[i].z += forceij.z
				force[j].z -= forceij.z
			}
		}
	}
	sim.pot = Epot
}
*/
