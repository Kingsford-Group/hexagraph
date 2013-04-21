import java.util.Vector;

// Our global data
int[] assignment;
int[] cellToNode;
Vector[] blockedSides;
String[] names;

Vector paths;
int n, gridx, gridy;
int hexRadius = 20;
float BLOCKED_OFFSET = 0.8;

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

int pos_mod(int i, int n) {
  return (i >= 0) ? (i % n) : (n + (i%n));
}

PVector[] drawHexSunken(float cx, float cy, float r, String t, Vector sides) {
  float a = 2*PI/6;
  float rp = BLOCKED_OFFSET * r;
  float s = r - rp;
  PVector[] out = new PVector[6];

  PVector[] outer = new PVector[6];
  PVector[] inner = new PVector[6];
  PVector[] iohalf = new PVector[12];
  for (int i = 0; i < 6; i++) {
    outer[i] = new PVector(r*cos(a*i) + cx, r*sin(a*i) + cy);
    inner[i] = new PVector(rp*cos(a*i) + cx, rp*sin(a*i) + cy);

    iohalf[pos_mod(2*i - 1, 12)] = new PVector(
    inner[i].x + s * cos(a*pos_mod(i-1, 6)), 
    inner[i].y + s * sin(a*pos_mod(i-1, 6))
      );
    iohalf[2*i] = new PVector(inner[i].x + s * cos(a*pos_mod(i+1,6)), inner[i].y + s * sin(a*pos_mod(i+1,6)));
  }

  String [] order = {
    "se", "s", "sw", "nw", "n", "ne"
  };
  
  println(sides);

  int lx=0,ly=0;
  // draw the shape
  beginShape();
  for (int i = 0; i < 6; i++) {
    String curside = order[pos_mod(i, 6)];
    String prevside = order[pos_mod(i-1, 6)];

    print(i);
    if (sides.contains(curside) && sides.contains(prevside)) {
      vertex(inner[i%6].x, inner[i%6].y);
      print("I");
    } else if(sides.contains(curside)) {
      PVector P = iohalf[pos_mod(2 * i - 1, 12)];
      vertex(P.x, P.y);
      print("A");
    } else if(sides.contains(prevside)) {
      PVector P = iohalf[pos_mod(2*i, 12)];
      vertex(P.x, P.y);
      print("B");
    } else {
      vertex(outer[i%6].x, outer[i%6].y);
      print("O");
    }
  }
  endShape(CLOSE);

  // draw label
  pushStyle();
  fill(0);
  textAlign(CENTER);
  text(t, cx, cy+textAscent()/2);
  popStyle();

  return new PVector[0];
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
  assert(abs(sidelen - rad) < 1e-4);
  print("C="); 
  println(sidelen);

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

  // if side is blocked, move BLOCKED_OFFSET pixels toward center
  if (cellToNode[c] >= 0) {
    println(c);
    Vector s = blockedSides[cellToNode[c]];
    for (int j = 0; j < s.size(); j++) {
      if (side.equals((String)s.get(j))) {
        xy[0] = (int)(hexCenters[c][0] + (xy[0] - hexCenters[c][0]) * BLOCKED_OFFSET);
        xy[1] = (int)(hexCenters[c][1] + (xy[1] - hexCenters[c][1]) * BLOCKED_OFFSET);
        break;
      }
    }
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
      assert(s[1].equals("hex"));
      gridx = int(s[2]);
      gridy = int(s[3]);
      n = int(s[4]);

      cellToNode = new int[gridx * gridy];
      for (int j = 0; j < cellToNode.length; j++) cellToNode[j] = -1;

      assignment = new int[n];
      names = new String[n];
      paths = new Vector();
      blockedSides = new Vector[n];
      break;

    case 'A':
      int u = int(s[1]);
      assignment[u] = int(s[2]);
      cellToNode[int(s[2])] = u;
      names[u] = s[3];

      // read the blocked sides into a vector
      blockedSides[u] = new Vector();
      String[] p = split(s[4], ',');
      for (int j = 0; j < p.length; j++) blockedSides[u].add(p[j]);
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
    drawHexSunken(cx, cy, hexRadius, names[i], blockedSides[i]); //+ "," + str(assignment[i]));
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
    if (C.substring(0, 2).equals("c=")) {
      String[] c = split(C.substring(2), ",");
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
  size(1000, 500);
  PFont f = createFont("Helvetica", 16, true);
  textFont(f, 12);

  readHexLayout("hex.layout");
  int[][] hexCenters = createHexCenters(gridx, gridy, 20);

  drawNodes(assignment, hexCenters);
  drawPaths(paths, hexCenters);
  noFill();
}

