#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>

#include "Garfield/ComponentElmer2d.hh"

namespace {

void PrintErrorReadingFile(const std::string& hdr, const std::string& file,
                           const int line) {
  std::cerr << hdr << "\n    Error reading file " << file << " (line " << line
            << ").\n";
}

}

namespace Garfield {

ComponentElmer2d::ComponentElmer2d() : ComponentFieldMap("Elmer2d") {
  m_is3d = false;
  m_elementType = ElementType::Serendipity;

  // Default bounding box
  m_minBoundingBox[2] = -50;
  m_maxBoundingBox[2] = 50;
}

ComponentElmer2d::ComponentElmer2d(const std::string& header,
                                   const std::string& elist,
                                   const std::string& nlist,
                                   const std::string& mplist,
                                   const std::string& volt, 
                                   const std::string& unit)
    : ComponentElmer2d() {

  Initialise(header, elist, nlist, mplist, volt, unit);
}

bool ComponentElmer2d::Initialise(const std::string& header,
                                const std::string& elist,
                                const std::string& nlist,
                                const std::string& mplist,
                                const std::string& volt,
                                const std::string& unit) {
  const std::string hdr = m_className + "::Initialise:";
  Reset();

  // Keep track of the success.
  bool ok = true;

  // Buffer for reading
  constexpr int size = 100;
  char line[size];

  // Open the header.
  std::ifstream fheader(header);
  if (!fheader) {
    PrintCouldNotOpen("Initialise", header);
    return false;
  }

  // Temporary variables for use in file reading
  char* token = nullptr;
  bool readerror = false;
  bool readstop = false;
  int il = 0;

  // Read the header to get the number of nodes and elements.
  fheader.getline(line, size, '\n');
  token = strtok(line, " ");
  const int nNodes = ReadInteger(token, 0, readerror);
  token = strtok(nullptr, " ");
  const int nElements = ReadInteger(token, 0, readerror);
  std::cout << hdr << "\n    Read " << nNodes << " nodes and " << nElements
            << " elements from file " << header << ".\n";
  if (readerror) {
    PrintErrorReadingFile(hdr, header, il);
    fheader.close();
    return false;
  }

  // Close the header file.
  fheader.close();

  // Open the nodes list.
  std::ifstream fnodes(nlist);
  if (!fnodes) {
    PrintCouldNotOpen("Initialise", nlist);
    return false;
  }

  // Check the value of the unit.
  double funit = ScalingFactor(unit);
  if (funit <= 0.) {
    std::cerr << hdr << " Unknown length unit " << unit << " Will use cm.\n";
    funit = 1.0;
  }
  if (m_debug) std::cout << hdr << " Unit scaling factor = " << funit << ".\n";

  // Read the nodes from the file.
  for (il = 0; il < nNodes; il++) {
    // Get a line from the nodes file.
    fnodes.getline(line, size, '\n');

    // Ignore the first two characters.
    token = strtok(line, " ");
    token = strtok(nullptr, " ");

    // Get the node coordinates.
    token = strtok(nullptr, " ");
    double xnode = ReadDouble(token, -1, readerror);
    token = strtok(nullptr, " ");
    double ynode = ReadDouble(token, -1, readerror);
    token = strtok(nullptr, " ");
    double znode = ReadDouble(token, -1, readerror);
    if (readerror) {
      PrintErrorReadingFile(hdr, nlist, il);
      fnodes.close();
      return false;
    }

    // Set up and create a new node.
    Node node;
    node.w.clear();
    node.x = xnode * funit;
    node.y = ynode * funit;
    node.z = znode * funit;
    node.v = 0.;
    m_nodes.push_back(std::move(node));
  }

  // Close the nodes file.
  fnodes.close();

  // Open the potential file.
  std::ifstream fvolt(volt);
  if (!fvolt) {
    PrintCouldNotOpen("Initialise", volt);
    return false;
  }

  // Reset the line counter.
  il = 1;

  // Read past the header.
  while (!readstop && fvolt.getline(line, size, '\n')) {
    token = strtok(line, " ");
    if (strcmp(token, "Perm:") == 0) readstop = true;
    il++;
  }

  // Should have stopped: if not, print error message.
  if (!readstop) {
    std::cerr << hdr << "\n    Error reading past header of potentials file "
              << volt << ".\n";
    fvolt.close();
    return false;
  }

  // Read past the permutation map (number of lines = nNodes).
  for (int tl = 0; tl < nNodes; tl++) {
    fvolt.getline(line, size, '\n');
    il++;
  }

  // Read the potentials.
  for (int tl = 0; tl < nNodes; tl++) {
    fvolt.getline(line, size, '\n');
    token = strtok(line, " ");
    double v = ReadDouble(token, -1, readerror);
    if (readerror) {
      PrintErrorReadingFile(hdr, volt, il);
      fvolt.close();
      return false;
    }
    // Place the voltage in its appropriate node.
    m_nodes[tl].v = v;
  }

  // Close the potentials file.
  fvolt.close();

  // Open the materials file.
  std::ifstream fmplist(mplist);
  if (!fmplist) {
    PrintCouldNotOpen("Initialise", mplist);
    return false;
  }

  // Read the dielectric constants from the materials file.
  fmplist.getline(line, size, '\n');
  token = strtok(line, " ");
  if (readerror) {
    std::cerr << hdr << "\n    Error reading number of materials from "
              << mplist << ".\n";
    fmplist.close();
    return false;
  }
  const unsigned int nMaterials = ReadInteger(token, 0, readerror);
  m_materials.resize(nMaterials);
  for (auto& material : m_materials) {
    material.ohm = -1;
    material.eps = -1;
    material.medium = nullptr;
  }
  for (il = 2; il < ((int)nMaterials + 2); il++) {
    fmplist.getline(line, size, '\n');
    token = strtok(line, " ");
    ReadInteger(token, -1, readerror);
    token = strtok(nullptr, " ");
    double dc = ReadDouble(token, -1.0, readerror);
    if (readerror) {
      PrintErrorReadingFile(hdr, mplist, il);
      fmplist.close();
      return false;
    }
    m_materials[il - 2].eps = dc;
    std::cout << "    Set material " << il - 2 << " of "
              << nMaterials << " to eps " << dc << ".\n";
  }

  // Close the materials file.
  fmplist.close();

  // Find lowest epsilon, check for eps = 0, set default drift medium.
  if (!SetDefaultDriftMedium()) return false;

  // Open the elements file.
  std::ifstream felems(elist);
  if (!felems) {
    PrintCouldNotOpen("Initialise", elist);
    return false;
  }

  // Read the elements and their material indices.
  for (il = 0; il < nElements; il++) {
    // Get a line
    felems.getline(line, size, '\n');

    // Split into tokens.
    token = strtok(line, " ");
    // Read the 2nd-order element.
    token = strtok(nullptr, " ");
    int imat = ReadInteger(token, -1, readerror) - 1;
    token = strtok(nullptr, " ");
    std::vector<int> inode;
    for (size_t k = 0; k < 8; ++k) {
      token = strtok(nullptr, " ");
      const int in = ReadInteger(token, -1, readerror);
      if (!readerror) inode.push_back(in);
    }
    token = strtok(nullptr, " ");
    if (inode.size() != 8) {
      PrintErrorReadingFile(hdr, elist, il);
      std::cerr << "    Read " << inode.size() << " node indices for element"
                 << il << " (expected 8).\n";
      felems.close();
      return false;
    }

    if (m_debug && il < 10) {
      std::cout << "    Read nodes ";
      for (size_t k = 0; k < 8; ++k) {
        std::cout << inode[k];
        if (k < 7) std::cout << ", ";
      }
      std::cout << "\n    from element " << il + 1 << " of "
                << nElements << " with material " << imat << ".\n";
    }

    // Check the material number and ensure that epsilon is non-negative.
    if (imat < 0 || imat > (int)nMaterials) {
      std::cerr << hdr << "\n    Out-of-range material number on file " << elist
                << " (line " << il << ").\n";
      std::cerr << "    Element: " << il << ", material: " << imat << "\n";
      ok = false;
      break;
    }
    if (m_materials[imat].eps < 0) {
      std::cerr << hdr << "\n    Element " << il << " in " << elist << "\n"
                << "    uses material " << imat << " which does not have\n"
                << "    a positive permittivity in " << mplist << ".\n";
      ok = false;
      break;
    }

    // Check the node numbers.
    bool degenerate = false;
    for (size_t k = 0; k < 8; ++k) {
      if (inode[k] < 1) {
        std::cerr << hdr << "\n    Found a node number < 1 on file " << elist
                  << " (line " << il << ").\n    Element: " << il
                  << ", material: " << imat << ".\n";
        ok = false;
      } 
      // These elements must not be degenerate.
      for (size_t kk = k + 1; kk < 8; ++kk) {
        if (inode[k] == inode[kk]) degenerate = true;
      }
    }
    if (degenerate) {
      std::cerr << hdr << "\n    Element " << il << " of file " << elist
                << " is degenerate,\n"
                << "    no such elements are allowed in this type of map.\n";
      ok = false;
    }
    if (!ok) break;
    Element element;
    element.degenerate = false;
    // Store the material reference.
    element.matmap = imat;

    // Node references
    //               3 -- 6 -- 2
    //               |         |
    //               7         5
    //               |         |
    //               0 -- 4 -- 1
    // Note: Ordering of Elmer elements is described in the
    // ElmerSolver manual (appendix D. at the time of this comment).
    // In order to work properly with Coordinates4 and Coordinates5, the
    // Elmer ordering 0,1,2,3,4,5,6,7 (counter-clockwise about the element)
    // will need to be changed to 3,2,1,0,6,5,4,7 (clockwise about the element)   
    // if the normal of the defined element points in the -Z direction.

    // Check the direction of the element normal, +Z or -Z.
    double x01 = m_nodes[inode[1] - 1].x - m_nodes[inode[0] - 1].x;
    double y01 = m_nodes[inode[1] - 1].y - m_nodes[inode[0] - 1].y;
    double x12 = m_nodes[inode[2] - 1].x - m_nodes[inode[1] - 1].x;
    double y12 = m_nodes[inode[2] - 1].y - m_nodes[inode[1] - 1].y;
    const double crossprod = x01 * y12 - y01 * x12;
    if (crossprod < 0) {
      element.emap[3] = inode[0] - 1;
      element.emap[2] = inode[1] - 1;
      element.emap[1] = inode[2] - 1;
      element.emap[0] = inode[3] - 1;
      element.emap[6] = inode[4] - 1;
      element.emap[5] = inode[5] - 1;
      element.emap[4] = inode[6] - 1;
      element.emap[7] = inode[7] - 1;
    } else {
      element.emap[0] = inode[0] - 1;
      element.emap[1] = inode[1] - 1;
      element.emap[2] = inode[2] - 1;
      element.emap[3] = inode[3] - 1;
      element.emap[4] = inode[4] - 1;
      element.emap[5] = inode[5] - 1;
      element.emap[6] = inode[6] - 1;
      element.emap[7] = inode[7] - 1;
    }
    m_elements.push_back(std::move(element));
  }

  // Close the elements file.
  felems.close();
  if (!ok) return false;

  // Set the ready flag.
  m_ready = true;
  std::cout << hdr << " Finished.\n";

  // Remove weighting fields (if any).
  m_wfields.clear();
  m_wfieldsOk.clear();

  Prepare();
  return true;
}

bool ComponentElmer2d::SetWeightingField(const std::string& wvolt, 
                                         const std::string& label) {
  const std::string hdr = m_className + "::SetWeightingField:";
  if (!m_ready) {
    PrintNotReady("SetWeightingField");
    std::cerr << "    Weighting field cannot be added.\n";
    return false;
  }

  // Open the voltage list.
  std::ifstream fwvolt(wvolt);
  if (!fwvolt) {
    PrintCouldNotOpen("SetWeightingField", wvolt);
    return false;
  }

  // Check if a weighting field with the same label already exists.
  const size_t iw = GetOrCreateWeightingFieldIndex(label);
  if (iw + 1 != m_wfields.size()) {
    std::cout << m_className << "::SetWeightingField:\n"
              << "    Replacing existing weighting field " << label << ".\n";
  }
  m_wfieldsOk[iw] = false;

  // Temporary variables for use in file reading
  constexpr int size = 100;
  char line[size];
  char* token = nullptr;
  bool readerror = false;
  bool readstop = false;
  int il = 1;

  // Read past the header.
  while (!readstop && fwvolt.getline(line, size, '\n')) {
    token = strtok(line, " ");
    if (strcmp(token, "Perm:") == 0) readstop = true;
    il++;
  }

  // Should have stopped: if not, print error message.
  if (!readstop) {
    std::cerr << hdr << "\n    Error reading past header of potentials file "
              << wvolt << ".\n";
    fwvolt.close();
    return false;
  }

  // Read past the permutation map (number of lines = nNodes).
  const int nNodes = m_nodes.size();
  for (int tl = 0; tl < nNodes; tl++) {
    fwvolt.getline(line, size, '\n');
    il++;
  }

  // Read the potentials.
  for (int tl = 0; tl < nNodes; tl++) {
    double v;
    fwvolt.getline(line, size, '\n');
    token = strtok(line, " ");
    v = ReadDouble(token, -1, readerror);
    if (readerror) {
      PrintErrorReadingFile(hdr, wvolt, il);
      fwvolt.close();
      return false;
    }
    // Place the weighting potential at its appropriate node and index.
    m_nodes[tl].w[iw] = v;
  }

  // Close the potentials file.
  fwvolt.close();
  std::cout << hdr << "\n    Read potentials from file " << wvolt << ".\n";

  // Set the ready flag.
  m_wfieldsOk[iw] = true;
  return true;
}

void ComponentElmer2d::SetRangeZ(const double zmin, const double zmax) {
  if (fabs(zmax - zmin) <= 0.) {
    std::cerr << m_className << "::SetRangeZ: Zero range is not permitted.\n";
    return;
  }
  m_minBoundingBox[2] = m_mapmin[2] = std::min(zmin, zmax);
  m_maxBoundingBox[2] = m_mapmax[2] = std::max(zmin, zmax);
}

}  // namespace Garfield
