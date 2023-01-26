//
// ********************************************************************
// *                                                                  *
// * Copyright 2022 The TOPAS Collaboration                           *
// *                                                                  *
// * Permission is hereby granted, free of charge, to any person      *
// * obtaining a copy of this software and associated documentation   *
// * files (the "Software"), to deal in the Software without          *
// * restriction, including without limitation the rights to use,     *
// * copy, modify, merge, publish, distribute, sublicense, and/or     *
// * sell copies of the Software, and to permit persons to whom the   *
// * Software is furnished to do so, subject to the following         *
// * conditions:                                                      *
// *                                                                  *
// * The above copyright notice and this permission notice shall be   *
// * included in all copies or substantial portions of the Software.  *
// *                                                                  *
// * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,  *
// * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES  *
// * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND         *
// * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT      *
// * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,     *
// * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     *
// * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR    *
// * OTHER DEALINGS IN THE SOFTWARE.                                  *
// *                                                                  *
// ********************************************************************
//

// #include "TsParameterManager.hh"
#include "../parameter/TsParameterManager.hh"

#include "TsMagneticFieldMap.hh"
#include "TsVGeometryComponent.hh"

#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4Tokenizer.hh"
#include "G4ChordFinder.hh"

#include <fstream>
#include <locale>
#include <map>

// something something setting up the magnetic field
TsMagneticFieldMap::TsMagneticFieldMap(TsParameterManager* pM,TsGeometryManager* gM, TsVGeometryComponent* component):
TsVMagneticField(pM, gM, component), fInvertX(false), fInvertY(false), fInvertZ(false), fNX(0), fNY(0), fNZ(0) {
	ResolveParameters();
}

// maybe a function for cleaning everything up in a memory clearing situation?
// what does the ~ mean?
TsMagneticFieldMap::~TsMagneticFieldMap() {
	fFieldX.clear();
	fFieldY.clear();
	fFieldZ.clear();
	if(fChordFinder) delete fChordFinder;
}

// figure out the parameters of of the magnetic field we want
void TsMagneticFieldMap::ResolveParameters() {
	// open the file for reading in the values
	std::ifstream file(fPm->GetStringParameter(fComponent->GetFullParmName("MagneticField3DTable")));
	// if the file could not be opened (file has a value of zero)
	if (!file) {
		// output for being unable to open the given file
		G4cerr << "" << G4endl;
		G4cerr << "Topas is exiting due to a serious error." << G4endl;
		G4cerr << "The parameter: " << fComponent->GetFullParmName("MagneticField3DTable") << G4endl;
		G4cerr << "references a MagneticField3DTable file that cannot be found:" << G4endl;
		G4cerr << fPm->GetStringParameter(fComponent->GetFullParmName("MagneticField3DTable")) << G4endl;
		fPm->AbortSession(1);
	}
	// we now have an opened file with name file

	G4String line; // create a new string
	bool ReadingHeader = true;
	G4int counter = 0; // we are going to count up
	double xval = 0.,yval = 0.,zval = 0.,bx,by,bz; // create 6 doubles
	// an x, y, and z that are set to zero, then 3 more xyz values unintialized
	// these are the xval, yval, zval position information, corrosponding to the
	// bx, by, bz magnetic field information
	int ix = 0;
	int iy = 0;
	int iz = 0;
	// three integers to go with the three dimensions, to be used as indexs
	std::map<G4String,double> headerUnits; // something about comparisons with keys?
	// I think headerUnits somehow gets assigned the G4String and double things?
	// More recent Kepler's guess: This makes some sort of searchable binary tree
	// with items in it that have a value that is a double and a string, such as a
	// value and its units. You can search it by the string, and it returns the double?
	// code keeps dereferencing it by strings and getting a resulting number of some type,
	// probably a double
	std::vector<G4String> headerUnitStrings; // makes a headerUnitStrings vector.
	// basically a dynamically sized array that can have its size changed. Allows for
	// traveling to a location in the array by pointer offset, just as one would expect of
	// a typical array
	std::vector<G4String> headerFields; // same as above

	while (file.good()) {
		// while we still have lines of the file to read
		getline(file,line);
		// reads one line of the file and then puts it into "line"
		if (line.find_last_not_of(" \t\f\v\n\r") == std::string::npos)
			continue;
			// if the last non listed character is some sort of EOF that is given,
			// then we will leave this while loop

		std::string::size_type pos = line.find_last_not_of(' ');
		if(pos != std::string::npos) {
			// if pos is not this mysterious end of file operator
			line.erase(pos + 1);
			pos = line.find_first_not_of(" \t\n\f\v\r\n");
			if(pos != std::string::npos) line.erase(0, pos);
			// erase from 0 to pos?
		} else {
			line.erase(line.begin(), line.end());
			// otherwise erase the line from begin to end?
			// where are begin and end?
		}

		std::vector<G4String> thisRow;
		// vector of strings called thisRow

		G4Tokenizer next(line);
		// make something of class G4Tokenizer on line
		G4String token = next();
		// token is the result of next
		while (token != "" && token != "\t" && token != "\n" && token != "\r" && token != "\f" && token != "\v") {
			// while token is not the characters we don't care about
			thisRow.push_back(token);
			token = next();
			// push the letters onto the vector of thisRow
		}

		if ((thisRow[0] == "0") && (counter > 0)) {
			// Found end of header, signal start of data read
			if (headerUnitStrings.size() == 0) {
				if (headerFields.size() > 6) {
					// too many header fields, so more than 6 columns
					G4cerr << "" << G4endl;
					G4cerr << "Topas is exiting due to a serious error." << G4endl;
					G4cerr << "Header information was not usable from MagneticField3DTable file:" << G4endl;
					G4cerr << fPm->GetStringParameter(fComponent->GetFullParmName("MagneticField3DTable")) << G4endl;
					G4cerr << "Only six fields (x,y,z,Bx,By,Bz) are allowed without specified units. Please include explicit unit declaration in the header" << G4endl;
					fPm->AbortSession(1);
				} else {
					// we've hit the data and not found any header information
					// write default units to the header strings vector
					G4cout << "No units specified, setting to 'mm' for x,y,z and 'tesla' for Bx,By,Bz" << G4endl;
					headerUnitStrings.push_back("mm");
					headerUnitStrings.push_back("mm");
					headerUnitStrings.push_back("mm");
					headerUnitStrings.push_back("tesla");
					headerUnitStrings.push_back("tesla");
					headerUnitStrings.push_back("tesla");
				}
			}

			for(G4int i = 0; i < (G4int)headerFields.size(); i++) {
				// for i less than the number of headerFields
				G4String unitString = headerUnitStrings[i];
				// unitString has the value of the units of this column
				std::locale loc; // something to do with IO? 

				size_t f = unitString.find("[");
				// some sort of search result for G4Strings? Probably the index of the "[" char
				// there shouldn't be a "[" as we are past the header

				if (f != std::string::npos) {
					// if f has a usable value?
					unitString.replace(f, std::string("[").length(), "");
					// replace the opening bracket with a nothing
				};

				f = unitString.find("]");
				// find the closing bracket

				if (f != std::string::npos) {
					unitString.replace(f, std::string("]").length(), "");
					// replace the closing bracket with a nothing
				};

				for (std::string::size_type j = 0; j < unitString.length(); j++) {
					// for the length of the string
					unitString[j] = std::tolower(unitString[j],loc);
					// make everything lowercase
					// typecase friendliness
				}

				// set the units
				if (unitString == "mm") {
					headerUnits[headerFields[i]] = mm; // referencing by a string?
					// headerFields[i] should be a string!
					// and then it is set equal to mm? not "mm" or something sensable like that!
					// my guess is mm is a geant4 thing of some sort, a unit, or a
					// multiplier to get it into the common unit in the geant4 instance
				} else
					if (unitString == "m" || unitString == "metre" || unitString == "meter") {
						headerUnits[headerFields[i]] = m;
					} else
						if (unitString == "tesla") {
							headerUnits[headerFields[i]] = tesla;
						}
						else {
							headerUnits[headerFields[i]] = 1;
							// unitless?, or maybe just unit-conversion-less
						}
			}

			// set the field sizes to fNX, it must have been set by the header?
			fFieldX.resize(fNX);
			fFieldY.resize(fNX);
			fFieldZ.resize(fNX);
			// further resizing stuff
			for (int index_x=0; index_x<fNX; index_x++) {
				fFieldX[index_x].resize(fNY);
				fFieldY[index_x].resize(fNY);
				fFieldZ[index_x].resize(fNY);
				for (int index_y=0; index_y<fNY; index_y++) {
					fFieldX[index_x][index_y].resize(fNZ);
					fFieldY[index_x][index_y].resize(fNZ);
					fFieldZ[index_x][index_y].resize(fNZ);
				}
			}

			ReadingHeader = false;
			// we are now confident we are done with the header
			counter = 0;
			continue;
		}

		if (ReadingHeader) {
			if (counter == 0) {
				// read the numbers at the top of the file. Its probably the number
				// of values in x, y, and z based on the context
				fNX = atoi(thisRow[0]);
				fNY = atoi(thisRow[1]);
				fNZ = atoi(thisRow[2]);
			} else {
				if (thisRow.size() < 2) continue;
				// we have too little to read, skip

				headerFields.push_back(thisRow[1]);
				// add the second value in thisRow to the end of header fields
				if (thisRow.size() == 3)
					headerUnitStrings.push_back(thisRow[2]);
					// if the thisRow information has the expected size, add the
					// expected location of the header unit string to the
					// headerUnitStrings vector

				if (thisRow.size() > 3) {
					// the header has a larger than expected number of parts!
					G4cerr << "" << G4endl;
					G4cerr << "Topas is exiting due to a serious error." << G4endl;
					G4cerr << "Header information was not usable from MagneticField3DTable file:" << G4endl;
					G4cerr << fPm->GetStringParameter(fComponent->GetFullParmName("MagneticField3DTable")) << G4endl;
					G4cerr << "Header has an unknown format on line" << G4endl;
					G4cerr << line << G4endl;
					G4cerr << "This error can be triggered by mismatch of linux/windows end-of-line characters." << G4endl;
					G4cerr << "If the opera file was created in windows, try converting it with dos2unix" << G4endl;
					fPm->AbortSession(1);
				}
			}
		} else {
			// we are not reading the header!
			if (thisRow.size() != headerFields.size()) {
				// exactly what it says below
				G4cerr << "" << G4endl;
				G4cerr << "Topas is exiting due to a serious error." << G4endl;
				G4cerr << "Header information was not usable from MagneticField3DTable file:" << G4endl;
				G4cerr << fPm->GetStringParameter(fComponent->GetFullParmName("MagneticField3DTable")) << G4endl;
				G4cerr << "File contains columns not in the header." << G4endl;
				fPm->AbortSession(1);
			}

			// set the position values to the first three values in this line of the file!
			xval = atof(thisRow[0]);
			yval = atof(thisRow[1]);
			zval = atof(thisRow[2]);
			// set the B field values to the second three values in this line of the file!
			bx = atof(thisRow[3]);
			by = atof(thisRow[4]);
			bz = atof(thisRow[5]);

			if ( ix==0 && iy==0 && iz==0 ) {
				// if we are at index 0, aka this is our first time in this code segment
				// find the lowest value corner of the given set of position values
				fMinX = xval * headerUnits["X"];
				fMinY = yval * headerUnits["Y"];
				fMinZ = zval * headerUnits["Z"];
				// headerUnits is some sort of double containing a conversion factor
				// for each dim.
			}

			// calculate the field values in this location, and store them in the
			// field table!
			fFieldX[ix][iy][iz] = bx * headerUnits["BX"];
			fFieldY[ix][iy][iz] = by * headerUnits["BY"];
			fFieldZ[ix][iy][iz] = bz * headerUnits["BZ"];

			iz++;
			// move index up by 1 in the z direction
			if (iz == fNZ) {
				// if we have hit the maximum z value then move the y value by 1
				iy++;
				iz = 0;
			}

			if (iy == fNY) {
				// if we have hit the maximum y value then move the x value by 1
				ix++;
				iy = 0;
			}
		}

		counter++;
		// we have made it through one iteration of the main loading tick
	}

	// done with loading the file, we no longer need it open
	file.close();

	if (fNX == 0) {
		// a E field with zero values in x was given! what sillyness is this!
		G4cerr << "" << G4endl;
		G4cerr << "Topas is exiting due to a serious error." << G4endl;
		G4cerr << "Header information was not usable from MagneticField3DTable file:" << G4endl;
		G4cerr << fPm->GetStringParameter(fComponent->GetFullParmName("MagneticField3DTable")) << G4endl;
		fPm->AbortSession(1);
	}

	// we are now at the far corner of this box, so find the converted location
	fMaxX = xval * headerUnits["X"];
	fMaxY = yval * headerUnits["Y"];
	fMaxZ = zval * headerUnits["Z"];

	if (fMaxX < fMinX) {
		// if our values went from large to small in x
		std::swap(fMaxX,fMinX); // swap max and min corner values
		fInvertX = true; // set to inverted mode in x
	}

	if (fMaxY < fMinY) {
		// same as above but for y
		std::swap(fMaxY,fMinY);
		fInvertY = true;
	}

	if (fMaxZ < fMinZ) {
		// ibid but for z
		std::swap(fMaxZ,fMinZ);
		fInvertZ = true;
	}

	fDX = fMaxX - fMinX;
	fDY = fMaxY - fMinY;
	fDZ = fMaxZ - fMinZ;
	// total distance across the box in x, y, and z

	const G4RotationMatrix* rotM = fComponent->GetRotRelToWorld();
	// define a rotation matrix
	G4Point3D* fTransRelToWorld = GetComponent()->GetTransRelToWorld();
	// translate based on the translation matrix?
	G4ThreeVector transl = G4ThreeVector(fTransRelToWorld->x(),fTransRelToWorld->y(),fTransRelToWorld->z());
	// vector for translation!
	fAffineTransf = G4AffineTransform(rotM,transl);
	// some sort of transform
}


// now the function that actually gets called by geant4 to get the field
void TsMagneticFieldMap::GetFieldValue(const G4double Point[3], G4double* Field) const {
	const G4ThreeVector localPoint = fAffineTransf.Inverse().TransformPoint(G4ThreeVector(Point[0],Point[1],Point[2]));
	// my guess is that localPoint is the actual location we are being asked about

	// Tabulated 3D table has it's own field area and the region is supposed to be smaller than volume
	// Therefore it is necessary to check that the point is inside the field area.
	G4double FieldX;
	G4double FieldY;
	G4double FieldZ;

	if ( localPoint.x() >= fMinX && localPoint.x() <= fMaxX && localPoint.y() >= fMinY && localPoint.y() <=fMaxY && localPoint.z() >=fMinZ && localPoint.z() <= fMaxZ ) {
		// Position of given point within region, normalized to the range [0,1]
		G4double xFraction = (localPoint.x() - fMinX)/fDX;
		G4double yFraction = (localPoint.y() - fMinY)/fDY;
		G4double zFraction = (localPoint.z() - fMinZ)/fDZ;

		if (fInvertX)
			xFraction = 1 - xFraction;
		if (fInvertY)
			yFraction = 1 - yFraction;
		if (fInvertZ)
			zFraction = 1 - zFraction;

		// Position of the point within the cuboid defined by the
		// nearest surrounding tabulated points
		G4double xDIndex; // future locations of the indexs
		G4double yDIndex;
		G4double zDIndex;
		// some sort of mod function giving a pair of doubles. The returned value
		// xLocal is the part after the decimal place, xDIndex is the integer part
		G4double xLocal = ( std::modf(xFraction*(fNX-1), &xDIndex));
		G4double yLocal = ( std::modf(yFraction*(fNY-1), &yDIndex));
		G4double zLocal = ( std::modf(zFraction*(fNZ-1), &zDIndex));

		// The indices of the nearest tabulated point whose coordinates
		// are all less than those of the given point
		G4double xIndex = static_cast<int>(xDIndex);
		G4double yIndex = static_cast<int>(yDIndex);
		G4double zIndex = static_cast<int>(zDIndex);

		// In rare cases, value is all the way to the end of the first bin.
		// Need to make sure it is assigned to that bin and not to the non-existant next bin.
//		if (xIndex + 1 == fNX) xIndex--;
		if (xIndex + 1 == fNX) {
			xIndex--;
			xLocal = 1;
		}
//		if (yIndex + 1 == fNY) yIndex--;
		if (yIndex + 1 == fNY) {
			yIndex--;
			yLocal = 1;
		}
//		if (zIndex + 1 == fNZ) zIndex--;
		if (zIndex + 1 == fNZ) {
			zIndex--;
			zLocal = 1;
		}

		// Full 3-dimensional version
		FieldX =
		fFieldX[xIndex  ][yIndex  ][zIndex  ] * (1-xLocal) * (1-yLocal) * (1-zLocal) +
		fFieldX[xIndex  ][yIndex  ][zIndex+1] * (1-xLocal) * (1-yLocal) *    zLocal  +
		fFieldX[xIndex  ][yIndex+1][zIndex  ] * (1-xLocal) *    yLocal  * (1-zLocal) +
		fFieldX[xIndex  ][yIndex+1][zIndex+1] * (1-xLocal) *    yLocal  *    zLocal  +
		fFieldX[xIndex+1][yIndex  ][zIndex  ] *    xLocal  * (1-yLocal) * (1-zLocal) +
		fFieldX[xIndex+1][yIndex  ][zIndex+1] *    xLocal  * (1-yLocal) *    zLocal  +
		fFieldX[xIndex+1][yIndex+1][zIndex  ] *    xLocal  *    yLocal  * (1-zLocal) +
		fFieldX[xIndex+1][yIndex+1][zIndex+1] *    xLocal  *    yLocal  *    zLocal;
		FieldY =
		fFieldY[xIndex  ][yIndex  ][zIndex  ] * (1-xLocal) * (1-yLocal) * (1-zLocal) +
		fFieldY[xIndex  ][yIndex  ][zIndex+1] * (1-xLocal) * (1-yLocal) *    zLocal  +
		fFieldY[xIndex  ][yIndex+1][zIndex  ] * (1-xLocal) *    yLocal  * (1-zLocal) +
		fFieldY[xIndex  ][yIndex+1][zIndex+1] * (1-xLocal) *    yLocal  *    zLocal  +
		fFieldY[xIndex+1][yIndex  ][zIndex  ] *    xLocal  * (1-yLocal) * (1-zLocal) +
		fFieldY[xIndex+1][yIndex  ][zIndex+1] *    xLocal  * (1-yLocal) *    zLocal  +
		fFieldY[xIndex+1][yIndex+1][zIndex  ] *    xLocal  *    yLocal  * (1-zLocal) +
		fFieldY[xIndex+1][yIndex+1][zIndex+1] *    xLocal  *    yLocal  *    zLocal;
		FieldZ =
		fFieldZ[xIndex  ][yIndex  ][zIndex  ] * (1-xLocal) * (1-yLocal) * (1-zLocal) +
		fFieldZ[xIndex  ][yIndex  ][zIndex+1] * (1-xLocal) * (1-yLocal) *    zLocal  +
		fFieldZ[xIndex  ][yIndex+1][zIndex  ] * (1-xLocal) *    yLocal  * (1-zLocal) +
		fFieldZ[xIndex  ][yIndex+1][zIndex+1] * (1-xLocal) *    yLocal  *    zLocal  +
		fFieldZ[xIndex+1][yIndex  ][zIndex  ] *    xLocal  * (1-yLocal) * (1-zLocal) +
		fFieldZ[xIndex+1][yIndex  ][zIndex+1] *    xLocal  * (1-yLocal) *    zLocal  +
		fFieldZ[xIndex+1][yIndex+1][zIndex  ] *    xLocal  *    yLocal  * (1-zLocal) +
		fFieldZ[xIndex+1][yIndex+1][zIndex+1] *    xLocal  *    yLocal  *    zLocal;
		// these give the values of the field at point x.,y,z, doing a linear mixing
		// of the values of the 8 known field values at the points around where we are
		// looking

		G4ThreeVector B_local = G4ThreeVector(FieldX,FieldY,FieldZ);
		G4ThreeVector B_global = fAffineTransf.TransformAxis(B_local);
		// make a B field vector in local space, then transform it into global space

		Field[0] = B_global.x() ;
		Field[1] = B_global.y() ;
		Field[2] = B_global.z() ;
		// assign the field values
	} else {
		Field[0] = 0.0;
		Field[1] = 0.0;
		Field[2] = 0.0;
		// give zero field from this outside if it was outside of the box we know
	}
}
