#include <math.h>
#include <string>
#include <assert.h>

#include "lipid.h"
#include "membrane.h"
#include "periodicboundary.h"


//
//     DPPC  DPPC  DIPC  DPPC     DPPC  DPPC  DIPC  DPPC     DPPC  DPPC  DIPC  DPPC
//  DIPC  DPPC  DIPC  DPPC     DIPC  DPPC  DIPC  DPPC     DIPC  DPPC  DIPC  DPPC
//                            |-------------------------|
//     DPPC  DPPC  DIPC  DPPC |   DPPC  DPPC  DIPC  DPPC|    DPPC  DPPC  DIPC  DPPC
//  DIPC  DPPC  DIPC  DPPC    |DIPC  DPPC  DIPC  DPPC   | DIPC  DPPC  DIPC  DPPC
//                            |-------------------------|
//     DPPC  DPPC  DIPC  DPPC     DPPC  DPPC  DIPC  DPPC     DPPC  DPPC  DIPC  DPPC
//  DIPC  DPPC  DIPC  DPPC     DIPC  DPPC  DIPC  DPPC     DIPC  DPPC  DIPC  DPPC
//


// returns pbc-corrected grid coordinates of neighbours on the hex-grid
// 0: top-left, 1: top-right, 2: left, 3: right, 4:bottom-left, 5: bottom-right
void periodic_neighbours(membrane& leaflet, int i, int j, int nbs[][2]) {
	
	int n = leaflet.getgrid().size();

	nbs[2][0] = i;   nbs[2][1] = j-1;
	nbs[3][0] = i;   nbs[3][1] = j+1;
	if (i % 2 == 0) {
		nbs[0][0] = i - 1; nbs[0][1] = j;
		nbs[1][0] = i - 1; nbs[1][1] = j + 1;
		nbs[4][0] = i + 1; nbs[4][1] = j;
		nbs[5][0] = i + 1; nbs[5][1] = j + 1;
	}
	else {
		nbs[0][0] = i - 1; nbs[0][1] = j - 1;
		nbs[1][0] = i - 1; nbs[1][1] = j;
		nbs[4][0] = i + 1; nbs[4][1] = j - 1;
		nbs[5][0] = i + 1; nbs[5][1] = j;
	}
	for (int k = 0; k < 6; k++) {
		if (nbs[k][0] < 0) {
			nbs[k][0] = n - 1;
		}
		if (nbs[k][0] == n) {
			nbs[k][0] = 0;
		}
		if (nbs[k][1] < 0) {
			nbs[k][1] = n - 1;
		}
		if (nbs[k][1] == n) {
			nbs[k][1] = 0;
		}
	}
}