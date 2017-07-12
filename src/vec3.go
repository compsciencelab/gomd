
type Real float32

/*
type Vec3 struct  {
	x,y,z Real
}
*/

type Vec3 [3]Real


func add(a,b Vec3) Vec3 {
	return Vec3(a[0]+b[0], a[1]+b[1],a[2]+b[2]])
}