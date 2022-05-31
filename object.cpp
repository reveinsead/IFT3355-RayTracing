#include "object.hpp"

#include <cmath>
#include <cfloat>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <iostream>
#include <algorithm> 


bool Object::intersect(Ray ray, Intersection& hit) const
{
	// Assure une valeur correcte pour la coordonnée W de l'origine et de la direction
	// Vous pouvez commentez ces lignes si vous faites très attention à la façon de construire vos rayons.
	ray.origin[3] = 1;
	ray.direction[3] = 0;

	Ray local_ray(i_transform * ray.origin, i_transform * ray.direction);
	//!!! NOTE UTILE : pour calculer la profondeur dans localIntersect(), si l'intersection se passe à
	// ray.origin + ray.direction * t, alors t est la profondeur
	//!!! NOTE UTILE : ici, la direction peut êytre mise à l'échelle, alors vous devez la renormaliser
	// dans localIntersect(), ou vous aurez une profondeur dans le système de coordonnées local, qui
	// ne pourra pas être comparée aux intersection avec les autres objets.
	if (localIntersect(local_ray, hit))
	{
		// Assure la valeur correcte de W.
		hit.position[3] = 1;
		hit.normal[3] = 0;

		// Transforme les coordonnées de l'intersection dans le repère global.
		hit.position = transform * hit.position;
		hit.normal = (n_transform * hit.normal).normalized();

		return true;
	}

	return false;
}


bool Sphere::localIntersect(Ray const& ray, Intersection& hit) const
{
	// @@@@@@ VOTRE CODE ICI
	// Vous pourriez aussi utiliser des relations géométriques pures plutôt que les
	// outils analytiques présentés dans les slides.
	// Ici, dans le système de coordonées local, la sphère est centrée en (0, 0, 0)
	//
	// NOTE : hit.depth est la profondeur de l'intersection actuellement la plus proche,
	// donc n'acceptez pas les intersections qui occurent plus loin que cette valeur.
	double t;
	double A = ray.direction.dot(ray.direction);
	double B = ray.origin.dot(ray.direction) * 2;
	double C = ray.origin.dot(ray.origin) - pow(this->radius, 2);
	double D = pow(B, 2) - 4 * A * C;//determinant

	if (D < 0) {
		return false;
	}
	t = (-B - sqrt(D)) / (2 * A);//always take the min root

	if (t < 0) {
		return false;
	}
	hit.position = ray.origin + t * ray.direction;
	hit.depth = t;
	//update values
	Vector outWardNormal = hit.position / this->radius;
	if (ray.direction.dot(outWardNormal) > 0) {
		hit.normal = -outWardNormal;
		return true;
	}
	else {
		hit.normal = outWardNormal;
		return true;
	}

}

bool Plane::localIntersect(Ray const& ray, Intersection& hit) const
{
	// @@@@@@ VOTRE CODE ICI
	// N'acceptez pas les intersections tant que le rayon est à l'intérieur du plan.
	// ici, dans le système de coordonées local, le plan est à z = 0.
	//
	// NOTE : hit.depth est la profondeur de l'intersection actuellement la plus proche,
	// donc n'acceptez pas les intersections qui occurent plus loin que cette valeur.


	//Method: P0=any point on the plane,p=hit.position
	//(p - p0)*n = 0; Plane's normal is perpendicular to any vector on the plane
	//p=ray.origin+ray.direction*t;
	//(ray.origin+ray.direction*t - p0)*n = 0
	//t=(p0-ray.origin)*n/ray.direction*n
	Vector planeNormal = Vector(0, 0, 1);
	Vector p0 = Vector(0, 0, 0);
	double t = (p0 - ray.origin).dot(planeNormal) / ray.direction.dot(planeNormal);

	//special case when line and plane are parallel but not contained in plane
	if (ABS_FLOAT(ray.direction.dot(planeNormal)) < 1e-6) {
		if (!(ABS_FLOAT((p0 - ray.origin).dot(planeNormal)) < 1e-6)) {
			return false;//no intersection
		}
		else {
			t = 0;//intersection along the line
		}
	}
	//if ray.direction*normal==0,ray is parallel to the plane,no intersection
	if (t < 0 || t>hit.depth) {
		return false;
	}

	Vector hitPoint = ray.origin + ray.direction * t;
	//update values
	hit.position = hitPoint;
	hit.depth = t;

	Vector outWardNormal = planeNormal;
	if (ray.direction.dot(outWardNormal) > 0) {
		hit.normal = -outWardNormal;
		return true;
	}
	else {
		hit.normal = outWardNormal;
		return true;
	}
}


bool Conic::localIntersect(Ray const& ray, Intersection& hit) const {
	// @@@@@@ VOTRE CODE ICI (licence créative)

	double zMin = this->zMin;
	double zMax = this->zMax;
	double radius1 = this->radius1;
	double radius2 = this->radius2;
	double A, B, C, D;
	double t;

	//if one of the radius is 0,it's a unit cone
	if (radius1 == 0) {
		radius1 = radius2;
	}
	if (radius2 == 0) {
		radius2 = radius1;
	}
	//case non unit cone:x^2/a^2+u^2/b^2=z^2
	//plugin x=x0+t*xDir,y=y0+t*yDir,z=z0+t*zDir
	A = pow(ray.direction[0], 2) / pow(radius1, 2) + pow(ray.direction[1], 2) / pow(radius2, 2) - pow(ray.direction[2], 2);
	B = 2 * ray.origin[0] * ray.direction[0] / pow(radius1, 2) + 2 * ray.origin[1] * ray.direction[1] / pow(radius2, 2) - 2 * ray.origin[2] * ray.direction[2];
	C = pow(ray.origin[0], 2) / pow(radius1, 2) + pow(ray.origin[1], 2) / pow(radius2, 2) - pow(ray.origin[2], 2);
	D = pow(B, 2) - 4 * A * C;//determinant

	//only one intersection point
	if (ABS_FLOAT(A) < 1e-6) {
		t = -C / B;
	}
	//no intersection point
	if (D < 0) {
		return false;
	}

	double t0 = (-B - sqrt(D)) / (2 * A);//compute the min root
	double t1 = (-B + sqrt(D)) / (2 * A);//compute the max root
	double z0 = ray.origin[2] + t0 * ray.direction[2];
	double z1 = ray.origin[2] + t1 * ray.direction[2];
	//check if z0 and z1 meets the lower bound and upperbound
	if ((z0 > zMax || z0 < zMin) && (z1 > zMax || z1 < zMin)) {
		return false;
	}
	else if ((z0 > zMax || z0 < zMin)) {
		t = t1;//case z0 is not in the interval[zMin,zMax]
	}
	else if ((z1 > zMax || z1 < zMin)) {
		t = t0;//case z1 is not in the interval[zMin,zMax]
	}
	else {
		t = std::min(t0, t1);//case two solutions
	}

	if (t < 0) {
		return false;
	}

	hit.position = ray.origin + t * ray.direction;
	hit.depth = t;
	//partial derivative of the implicite function
	Vector outWardNormal = Vector(2 * hit.position[0] / pow(radius1, 2), 2 * hit.position[1] / pow(radius2, 2), -2 * hit.position[2]).normalized();
	if (ray.direction.dot(outWardNormal) > 0) {
		hit.normal = -outWardNormal;
		return true;
	}
	else {
		hit.normal = outWardNormal;
		return true;
	}


}


// Intersections !
bool Mesh::localIntersect(Ray const& ray, Intersection& hit) const
{
	// Test de la boite englobante
	double tNear = -DBL_MAX, tFar = DBL_MAX;
	for (int i = 0; i < 3; i++) {
		if (ray.direction[i] == 0.0) {
			if (ray.origin[i] < bboxMin[i] || ray.origin[i] > bboxMax[i]) {
				// Rayon parallèle à un plan de la boite englobante et en dehors de la boite
				return false;
			}
			// Rayon parallèle à un plan de la boite et dans la boite: on continue
		}
		else {
			double t1 = (bboxMin[i] - ray.origin[i]) / ray.direction[i];
			double t2 = (bboxMax[i] - ray.origin[i]) / ray.direction[i];
			if (t1 > t2) std::swap(t1, t2); // Assure t1 <= t2

			if (t1 > tNear) tNear = t1; // On veut le plus lointain tNear.
			if (t2 < tFar) tFar = t2; // On veut le plus proche tFar.

			if (tNear > tFar) return false; // Le rayon rate la boite englobante.
			if (tFar < 0) return false; // La boite englobante est derrière le rayon.
		}
	}
	// Si on arrive jusqu'ici, c'est que le rayon a intersecté la boite englobante.

	// Le rayon interesecte la boite englobante, donc on teste chaque triangle.
	bool isHit = false;
	for (size_t tri_i = 0; tri_i < triangles.size(); tri_i++) {
		Triangle const& tri = triangles[tri_i];

		if (intersectTriangle(ray, tri, hit)) {
			isHit = true;
		}
	}
	return isHit;
}


double Mesh::implicitLineEquation(double p_x, double p_y,
	double e1_x, double e1_y,
	double e2_x, double e2_y) const
{
	return (e2_y - e1_y) * (p_x - e1_x) - (e2_x - e1_x) * (p_y - e1_y);
}

bool Mesh::intersectTriangle(Ray const& ray,
	Triangle const& tri,
	Intersection& hit) const
{
	// Extrait chaque position de sommet des données du maillage.
	Vector const& p0 = positions[tri[0].pi];
	Vector const& p1 = positions[tri[1].pi];
	Vector const& p2 = positions[tri[2].pi];

	// @@@@@@ VOTRE CODE ICI
	// Décidez si le rayon intersecte le triangle (p0,p1,p2).
	// Si c'est le cas, remplissez la structure hit avec les informations
	// de l'intersection et renvoyez true.
	// Vous pourriez trouver utile d'utiliser la routine implicitLineEquation()
	// pour calculer le résultat de l'équation de ligne implicite en 2D.
	//
	// NOTE : hit.depth est la profondeur de l'intersection actuellement la plus proche,
	// donc n'acceptez pas les intersections qui occurent plus loin que cette valeur.
	//!!! NOTE UTILE : pour le point d'intersection, sa normale doit satisfaire hit.normal.dot(ray.direction) < 0

	Vector p01 = p1 - p0;
	Vector p02 = p2 - p0;
	Vector triNormal = p01.cross(p02).normalized();

	double t = (p0 - ray.origin).dot(triNormal) / ray.direction.dot(triNormal);//see method in plane intersect
	//special case when line and plane are parallel but not contained in plane
	if (ABS_FLOAT(ray.direction.dot(triNormal)) < 1e-6) {
		if (!(ABS_FLOAT((p0 - ray.origin).dot(triNormal)) < 1e-6)) {
			return false;//no intersection
		}
		else {
			t = 0;//intersection along the line
		}
	}
	if (t > hit.depth || t < 0) {
		return false;
	}


	Vector hitPoint = ray.origin + ray.direction * t;

	//verify if the hitpoint is inside of a  triangle
	//fist calculate vectors 0->1,1->2,2->0
	Vector p12 = p2 - p1;
	Vector p20 = p0 - p2;
	//second calculate vectors 0->hitpoint,1->hitpoint,and 2->hitpoint
	Vector p0Hitpoint = hitPoint - p0;
	Vector p1Hitpoint = hitPoint - p1;
	Vector p2Hitpoint = hitPoint - p2;
    

	double area = p01.cross(p02).length();
	double alpha = p01.cross(p0Hitpoint).length() / area;
	double beta = p12.cross(p1Hitpoint).length() / area;
	double gamma = p20.cross(p2Hitpoint).length() / area;

	if ((alpha + gamma + beta - 1) > 1e-6) {
		return false;
	}
	

	//update values
	hit.position = hitPoint;
	hit.depth = t;
	hit.normal = triNormal;
	return true;
	Vector outWardNormal = triNormal;
	if (ray.direction.dot(outWardNormal) > 0) {
		hit.normal = -outWardNormal;
		return true;
	}
	else {
		hit.normal = outWardNormal;
		return true;
	}


}
