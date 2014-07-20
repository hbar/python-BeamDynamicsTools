int pnpoly(int nvert, float *vertx, float *verty, float testx, float testy)
{
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return c;
}

Argument	Meaning
nvert 	Number of vertices in the polygon. Whether to repeat the first vertex at the end is discussed below.
vertx, verty 	Arrays containing the x- and y-coordinates of the polygon's vertices.
testx, testy	X- and y-coordinate of the test point.
