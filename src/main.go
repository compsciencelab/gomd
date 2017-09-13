package main

import (
	"fmt"
	"math"

	"./sim"
)

func main() {
	natoms := 3
	niter := 1000
	sim := sim.New(natoms, sim.Vec3{10, 10, 10})
	sim.Read("../input_coor.dat")
	sim.Write("test.dat")
	sim.ComputeNonBonded()
	fmt.Println(0, sim.Pot, sim.Kin, sim.Pot+sim.Kin)
	for n := 0; n < niter; n++ {
		sim.FirstVV()
		sim.ResetForce()
		sim.ComputeNonBonded()
		sim.SecondVV()
		temp := sim.GetTemperature()
		//sim.ThermostatBerendesen(vel,natoms,dt)
		if math.Remainder(float64(n), 1) == 0 {
			fmt.Println(n, sim.Pot, sim.Kin, sim.Kin+sim.Pot, temp)
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
