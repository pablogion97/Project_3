#include "Domain.h"
using namespace std;
constexpr double SMALL = 1e-5;

Domain::Domain(Curvebase& bottom, Curvebase& right, Curvebase& top, Curvebase& left) {

	// Copy pointers of the curves to the array of pointers
	sides[0] = &bottom;
	sides[1] = &right;
	sides[2] = &top;
	sides[3] = &left;

	// Check the consistency of the provided curves
	if (!check_consistency()) {
		sides[0] = sides[1] = sides[2] = sides[3] = nullptr;
		cerr << "No valid boundary curves provided." << endl;
		exit(EXIT_FAILURE);
	}

	nx = ny = 0;
	x_ = y_ = nullptr;
}

// Copy constructor
Domain::Domain(const Domain& other) {
	for (int i = 0; i < 4; i++)
		sides[i] = other.sides[i];

	nx = other.nx;
	ny = other.ny;
	
	// Use memcpy for quicker processing if the grid has already been generated
	if (other.x_ != nullptr && other.y_ != nullptr) {
		x_ = new double[(nx + 1) * (ny + 1)];
		y_ = new double[(nx + 1) * (ny + 1)];
		memcpy(x_, other.x_, (nx + 1) * (ny + 1) * sizeof(double));
		memcpy(y_, other.y_, (nx + 1) * (ny + 1) * sizeof(double));
	}
	else {
		nx = ny = 0;
		x_ = y_ = nullptr;
	}
}

// Destructor
Domain::~Domain() {
	// Free the dynamically allocated memory upon destruction
	if (nx > 0 && ny > 0) {
		delete[] x_;
		delete[] y_;
	}
}

Domain& Domain::operator=(const Domain& other) {
	// self-assignment check expected
	if (this != &other) {

		// Copy side pointers
		for (int i = 0; i < 4; i++)
			sides[i] = other.sides[i];

		// If the sizes don't match
		if (nx != other.nx || ny != other.ny) {
			delete[] x_;
			delete[] y_;

			nx = other.nx;
			ny = other.ny;

			// If the grid was previously initialized
			if (other.x_ != nullptr && other.y_ != nullptr) {
				x_ = new double[(nx + 1) * (ny + 1)];
				y_ = new double[(nx + 1) * (ny + 1)];
			}
			else {
				// If not, then set nullptr and 0 sizes
				nx = ny = 0;
				x_ = y_ = nullptr;
				return *this;
			}

			// Otherwise, copy the contents onto the current instance
			memcpy(x_, other.x_, (nx + 1) * (ny + 1) * sizeof(double));
			memcpy(y_, other.y_, (nx + 1) * (ny + 1) * sizeof(double));
		}
	}

	return *this;
}

// Check if the start and end points of the curves match
bool Domain::check_consistency() const {
	bool correct = true;

	// Lower left corner
	correct &= abs(sides[3]->x(0.0) - sides[0]->x(0.0)) < SMALL;
	correct &= abs(sides[3]->y(0.0) - sides[0]->y(0.0)) < SMALL;

	// Lower right corner
	correct &= abs(sides[0]->x(1.0) - sides[1]->x(0.0)) < SMALL;
	correct &= abs(sides[0]->y(1.0) - sides[1]->y(0.0)) < SMALL;

	// Upper right corner
	correct &= abs(sides[1]->x(1.0) - sides[2]->x(1.0)) < SMALL;
	correct &= abs(sides[1]->y(1.0) - sides[2]->y(1.0)) < SMALL;

	// Upper left corner
	correct &= abs(sides[2]->x(0.0) - sides[3]->x(1.0)) < SMALL;
	correct &= abs(sides[2]->y(0.0) - sides[3]->y(1.0)) < SMALL;

	return correct;
}


// Generate the grid
void Domain::generate_grid(const int& nX, const int& nY) {
	// Check on user input
	if (nX <= 0 || nY <= 0) {
		cerr << "ERROR: Invalid grid sizes provided." << endl;
		exit(EXIT_FAILURE);
	}

	// Delete existing grid
	if (nx > 0 && ny > 0) { // There exists already a grid!
		delete[] x_;
		delete[] y_;
	}

	// Set sizes and allocate memory
	nx = nX;
	ny = nY;
	x_ = new double[(nx + 1) * (ny + 1)];
	y_ = new double[(nx + 1) * (ny + 1)];

	// Grid spacings (per unit length)
	const double hx = 1.0 / nx;
	const double hy = 1.0 / ny;

	// xi and eta variables for better readability
	double xi, eta;

	// Arrays for storing the boundary points
	Points* botPts = new Points[nx + 1];
	Points* topPts = new Points[nx + 1];
	Points* leftPts = new Points[ny + 1];
	Points* rightPts = new Points[ny + 1];

	// Compute the points of the lower and upper boundaries
	for (int i = 0; i <= nx; i++) {
		xi = hx * i;
		botPts[i] = Points{ sides[0]->x(xi), sides[0]->y(xi) };
		topPts[i] = Points{ sides[2]->x(xi), sides[2]->y(xi) };
	}

	// Compute the points of the left and right boundaries
	for (int i = 0; i <= ny; i++) {
		eta = 1.0 + tanh(3 * (hy * i - 1.0)) / tanh(3);
		leftPts[i] = Points{ sides[3]->x(eta), sides[3]->y(eta) };
		rightPts[i] = Points{ sides[1]->x(eta), sides[1]->y(eta) };
	}

	
	// Helper functions (lambda expressions) for better readability
	auto phi0 = [](double s) {return 1.0 - s; };
	auto phi1 = [](double s) {return s; };

	//Transfinite grid interpolation
	for (int j = 0; j <= ny; j++) {
		for (int i = 0; i <= nx; i++) {
			xi = hx * i;
			eta = 1.0 + tanh(3*(hy * j - 1.0))/tanh(3);

			x_[j * (nx + 1) + i] =
				phi0(xi) * leftPts[j].x
				+ phi1(xi) * rightPts[j].x
				+ phi0(eta) * botPts[i].x
				+ phi1(eta) * topPts[i].x
				- phi0(xi) * phi0(eta) * botPts[0].x
				- phi1(xi) * phi0(eta) * rightPts[0].x
				- phi0(xi) * phi1(eta) * botPts[0].x
				- phi1(xi) * phi1(eta) * rightPts[0].x;

			y_[j * (nx + 1) + i] =
				  phi0(xi) * leftPts[j].y
				+ phi1(xi) * rightPts[j].y
				+ phi0(eta) * botPts[i].y
				+ phi1(eta) * topPts[i].y
				- phi0(xi) * phi0(eta) * leftPts[0].y
				- phi1(xi) * phi0(eta) * leftPts[0].y
				- phi0(xi) * phi1(eta) * topPts[0].y
				- phi1(xi) * phi1(eta) * topPts[0].y;
		}
	}

	// Free memory (not strictly necessary)
	delete[] botPts;
	delete[] rightPts;
	delete[] topPts;
	delete[] leftPts;
	return;
}


// Export the grid as a csv file
void Domain::print_grid() const {
	// Check if the grid has been generated
	if (x_ == nullptr || y_ == nullptr) {
		cerr << "ERROR: No grid has been generated previously." << endl;
		exit(EXIT_FAILURE);
	}

	// Write grid to file as a csv table
	ofstream gridFile;
	gridFile.open("grid.csv");

	for (int j = 0; j <= ny; j++)
		for (int i = 0; i <= nx; i++)
			gridFile << x_[j * (nx + 1) + i] << "," << y_[j * (nx + 1) + i] << endl;

	gridFile.close();
}