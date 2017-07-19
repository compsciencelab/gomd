package main

import (
	"fmt"
	"math"

	"./sim"
)

func main() {
	natoms := 3188
	niter := 1000
	sim := sim.New(natoms, sim.Vec3{100, 100, 100})
	sim.Read("../input_coor.dat")
	sim.ComputeNonBonded()
	fmt.Println("ts Pot Kin Temp Tot")
	fmt.Println(0, sim.Pot, sim.Kin, 0, sim.Pot+sim.Kin)
	for n := 0; n < niter; n++ {
		sim.FirstVV()
		sim.ResetForce()
		sim.ComputeNonBonded()
		sim.SecondVV()
		temp := sim.GetTemperature()
		//sim.ThermostatBerendesen(vel,natoms,dt)
		if math.Remainder(float64(n), 10) == 0 {
			fmt.Println(n, sim.Kin, sim.Pot, temp, sim.Kin+sim.Pot)
		}
	}
	sim.Write("output.dat")
}
