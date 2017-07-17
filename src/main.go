package main

import (
	"fmt"
	"math"

	"./sim"
)

func main() {
	natoms := 10
	niter := 1
	sim := sim.New(natoms) //sim.New??
	//sim.ComputeNonBonded()
	//fmt.Println("%d %f %f %f\n", 0, Epot, Ekin, Epot+Ekin)
	for n := 0; n < niter; n++ {
		sim.FirstVV()
		//sim.force[:] = Vec3{0, 0, 0} //memset((void *)force, 0, natoms*sizeof(rvec));
		//sim.ComputeNonBonded()
		sim.SecondVV()
		temp := sim.GetTemperature()
		//sim.ThermostatBerendesen(vel,natoms,dt)
		if math.Remainder(float64(n), 10) == 0 {
			fmt.Println("%d %f %f %f %f %d\n", n, sim.Kin, sim.Pot, sim.Kin+sim.Pot, temp)
		}
	}
}

/*
func main() {
	natoms := 100
	//niter := 1
	sim := sim.New(natoms) //sim.New??
	fmt.Println("type ", sim.Pot)

}
*/
