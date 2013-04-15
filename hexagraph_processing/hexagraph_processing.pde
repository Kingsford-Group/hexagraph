import java.util.Vector;

// Our global data
int[] assignment;
String[] names;

Vector paths;
int n, gridx, gridy;
int hexRadius = 20;

//================================================================
// Hex Tiling Geometry
//================================================================

/* Draw a hexagon with the given radius "r" at center (cx,cy) */
void drawHex(float cx, float cy, float r, String t) {
  float a = 2*PI/6;
  beginShape();
  for (int i = 0; i <= 6; i++) {
    vertex(r*cos(a*i) + cx, r*sin(a*i) + cy);
  }
  endShape(CLOSE);
  pushStyle();
  fill(0);
  textAlign(CENTER);
  text(t, cx, cy+textAscent()/2);
  popStyle();
}

/* return the height of the hexagon from side to side */
float hexHeight(float hexRadius) {
  return 2*hexRadius * cos(2*PI / 12);
}

/* Construct a mapping from cell ids to hexgone centers */
int[][] createHexCenters(int gridw, int gridh, float rad) {
  int n = gridw * gridh; // number of cells
  int[][] hexCenters = new int[gridw*gridh][2]; // their locations

  // calculate geometry of a hexagon
  float hexheight = hexHeight(rad); 
  float sidelen = 2*rad * sin(2*PI / 12);

  // position every cell
  int x=0, y=0;
  int c = 0;
  for (int j = 0; j < gridh; ++j) {
    x = int((j % 2 == 1) ? (2*rad + sidelen / 2) : rad);
    y += int(hexheight / 2);

    for (int i = 0; i < gridw; ++i) {
      //drawHex(x,y,rad, str(c));
      hexCenters[c][0] = x;
      hexCenters[c][1] = y;
      c++;
      x += 3*rad;
    }
  }
  return hexCenters;
}


/* Get the xy coord of the center of a given side of cell c;
 side should be one of n,nw,ne,s,sw,se; must be lowercase */
int[] hexSidePosn(int c, String side, int[][] hexCenters) {
  int[] xy = new int[2];
  xy[0] = hexCenters[c][0];
  xy[1] = hexCenters[c][1];

  float hh = hexHeight(hexRadius) / 2;
  if (side.equals("s")) { 
    xy[1] += hh; 
    return xy;
  }
  if (side.equals("n")) { 
    xy[1] -= hh; 
    return xy;
  }

  float ydelta = cos(2 * PI / 6) * hh;
  float xdelta = sin(2 * PI / 6) * hh;

  if (side.charAt(0) == 'n') xy[1] -= ydelta;
  if (side.charAt(0) == 's') xy[1] += ydelta;
  if (side.length() >= 2) {
    if (side.charAt(1) == 'w') xy[0] -= xdelta;
    if (side.charAt(1) == 'e') xy[0] += xdelta;
  }

  return xy;
}

//================================================================
// Drawing a Hex Layout
//================================================================

/* Read the hex layout file */
void readHexLayout(String filename) {
  String[] lines = loadStrings(filename);
  for (int i = 0; i < lines.length; i++) {
    String[] s = split(lines[i], ' ');
    switch (s[0].charAt(0)) {
    case 'T':
      gridx = int(s[2]);
      gridy = int(s[3]);
      n = int(s[4]);

      assignment = new int[n];
      names = new String[n];
      paths = new Vector();
      break;

    case 'A':
      assignment[int(s[1])] = int(s[2]);
      names[int(s[1])] = s[3];
      break;

    case 'P':
      Vector a = new Vector();
      for (int j = 1; j < s.length; j++) a.add(s[j]);
      paths.add(a);
      break;
    }
  }
}


/* draw the cells that have nodes assigned to them */
void drawNodes(int[] assignment, int[][] hexCenters) {
  for (int i = 0; i < assignment.length; i++) {
    int cx = hexCenters[assignment[i]][0];
    int cy = hexCenters[assignment[i]][1];
    drawHex(cx, cy, hexRadius, names[i]); //+ "," + str(assignment[i]));
  }
}

/* draw curves along the given paths */
void drawPaths(Vector paths, int[][] hexCenters) {
  noFill();

  int[] xy = {
    0, 0
  };
  for (int i = 0; i < paths.size(); i++) {
    Vector P = (Vector)paths.get(i);
    assert(P.size() > 1);


    int start = 0;
    String C = (String)P.get(0);
    if (C.substring(0,2).equals("c=")) {
       String[] c = split(C.substring(2), ",");
       print(c[0]);
       assert(c.length == 3);
       stroke(float(c[0]), float(c[1]), float(c[2]));
       start = 1;
    }
    beginShape();
    for (int j = start; j < P.size(); j++) {
      String[] s = split((String)P.get(j), "_");
      assert(s.length == 2);
      xy = hexSidePosn(int(s[0]), s[1], hexCenters);
      curveVertex(xy[0], xy[1]);
      if (j == start) curveVertex(xy[0], xy[1]);
    }
    curveVertex(xy[0], xy[1]);
    endShape();
  }
}

void setup() {
  size(1000, 1000);
  PFont f = createFont("Helvetica", 16, true);
  textFont(f, 12);

  readHexLayout("hex.layout");
  int[][] hexCenters = createHexCenters(gridx, gridy, 20);

  drawNodes(assignment, hexCenters);
  drawPaths(paths, hexCenters);
  noFill();
}

