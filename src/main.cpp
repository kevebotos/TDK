#include <gmsh.h>
#include <Eigen/Sparse>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * Function to write temperature data to a Gmsh post-processing file (.pos)
 * @param filename: Name of the output file
 * @param node_coords: Vector containing node coordinates [x1, y1, z1, x2, y2, z2, ...]
 * @param temperatures: Vector containing temperature values for each node
 */
void write_pos_file(const std::string &filename,
                    const std::vector<double> &node_coords,
                    const std::vector<double> &temperatures)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << " for writing" << std::endl;
        return;
    }

    file << "View \"Temperature\" {\n";

    // Write node data as scalar points
    size_t num_nodes = temperatures.size();
    for (size_t i = 0; i < num_nodes; ++i)
    {
        file << "SP(" << node_coords[3 * i] << "," << node_coords[3 * i + 1] << "," << node_coords[3 * i + 2] << "){"
             << temperatures[i] << "};\n";
    }

    file << "};\n";
    file.close();

    std::cout << "Written temperature data to " << filename << std::endl;
}

/**
 * Function to calculate the area of a triangle given three vertices
 */
double triangle_area(double x1, double y1, double x2, double y2, double x3, double y3)
{
    return 0.5 * std::abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
}

/**
 * Function to compute local stiffness matrix for a triangular element
 * @param coords: Array of 6 doubles [x1, y1, x2, y2, x3, y3]
 * @param alpha: Thermal diffusivity
 * @return 3x3 local stiffness matrix
 */
Eigen::Matrix3d compute_local_stiffness(const double coords[6], double alpha)
{
    double x1 = coords[0], y1 = coords[1];
    double x2 = coords[2], y2 = coords[3];
    double x3 = coords[4], y3 = coords[5];

    // Calculate area
    double area = triangle_area(x1, y1, x2, y2, x3, y3);

    // Shape function gradients
    double b1 = y2 - y3, c1 = x3 - x2;
    double b2 = y3 - y1, c2 = x1 - x3;
    double b3 = y1 - y2, c3 = x2 - x1;

    Eigen::Matrix3d K_local;
    K_local.setZero();

    // Fill the stiffness matrix
    double factor = alpha / (4.0 * area);
    K_local(0, 0) = factor * (b1 * b1 + c1 * c1);
    K_local(0, 1) = factor * (b1 * b2 + c1 * c2);
    K_local(0, 2) = factor * (b1 * b3 + c1 * c3);
    K_local(1, 0) = K_local(0, 1);
    K_local(1, 1) = factor * (b2 * b2 + c2 * c2);
    K_local(1, 2) = factor * (b2 * b3 + c2 * c3);
    K_local(2, 0) = K_local(0, 2);
    K_local(2, 1) = K_local(1, 2);
    K_local(2, 2) = factor * (b3 * b3 + c3 * c3);

    return K_local;
}

/**
 * Function to determine material type based on element center
 * @param coords: Array of 6 doubles [x1, y1, x2, y2, x3, y3]
 * @return true if element is in burner region, false if in stove body
 */
bool is_burner_element(const double coords[6])
{
    // Calculate element center
    double cx = (coords[0] + coords[2] + coords[4]) / 3.0;
    double cy = (coords[1] + coords[3] + coords[5]) / 3.0;

    // Approximate burner locations (based on typical stove layout)
    // Burner centers: (3,1), (7,1), (3,5), (7,5) with radius ~0.5
    std::vector<std::pair<double, double>> burner_centers = {
        {3.0, 1.0}, // Burner 1
        {7.0, 1.0}, // Burner 2
        {3.0, 5.0}, // Burner 3
        {7.0, 5.0}  // Burner 4
    };

    double burner_radius = 0.6; // Radius to consider as burner region

    for (const auto &center : burner_centers)
    {
        double dx = cx - center.first;
        double dy = cy - center.second;
        double distance = std::sqrt(dx * dx + dy * dy);

        if (distance <= burner_radius)
        {
            return true;
        }
    }

    return false;
}

/**
 * Function to compute local mass matrix for a triangular element
 * @param coords: Array of 6 doubles [x1, y1, x2, y2, x3, y3]
 * @return 3x3 local mass matrix
 */
Eigen::Matrix3d compute_local_mass(const double coords[6])
{
    double x1 = coords[0], y1 = coords[1];
    double x2 = coords[2], y2 = coords[3];
    double x3 = coords[4], y3 = coords[5];

    // Calculate area
    double area = triangle_area(x1, y1, x2, y2, x3, y3);

    Eigen::Matrix3d M_local;
    M_local.setZero();

    // Mass matrix for linear triangular elements
    double diag_term = area / 6.0;      // 2/12 * area
    double off_diag_term = area / 12.0; // 1/12 * area

    M_local(0, 0) = diag_term;
    M_local(0, 1) = off_diag_term;
    M_local(0, 2) = off_diag_term;
    M_local(1, 0) = off_diag_term;
    M_local(1, 1) = diag_term;
    M_local(1, 2) = off_diag_term;
    M_local(2, 0) = off_diag_term;
    M_local(2, 1) = off_diag_term;
    M_local(2, 2) = diag_term;

    return M_local;
}

int main()
{
    // Simulation parameters
    double total_time = 300.0; // 5 minutes simulation
    double dt = 1.0;           // 1 second time step
    int output_interval = 10;  // Save a .pos file every 10 seconds

    // Material properties (realistic values)
    double alpha_stove = 4.2e-6;  // Thermal diffusivity for stainless steel (m²/s)
    double alpha_burner = 1.2e-5; // Thermal diffusivity for cast iron (m²/s)

    double initial_temp = 20.0;      // Room temperature (°C)
    double heat_source_temp = 200.0; // Target temperature for burner centers (°C)
    double ambient_temp = 20.0;      // Ambient temperature (°C)

    try
    {
        // Initialize Gmsh API
        gmsh::initialize();

        // Open the mesh file
        gmsh::open("stove.msh");

        std::cout << "Successfully opened stove.msh" << std::endl;

        // 1. Read node data
        std::vector<std::size_t> node_tags;
        std::vector<double> node_coords;
        std::vector<double> parametric_coords;

        gmsh::model::mesh::getNodes(node_tags, node_coords, parametric_coords);

        size_t num_nodes = node_tags.size();
        std::cout << "Number of nodes: " << num_nodes << std::endl;

        // Create a mapping from node tags to indices
        std::map<std::size_t, size_t> node_tag_to_index;
        for (size_t i = 0; i < num_nodes; ++i)
        {
            node_tag_to_index[node_tags[i]] = i;
        }

        // 2. Read triangular element data
        std::vector<int> element_types;
        std::vector<std::vector<std::size_t>> element_tags;
        std::vector<std::vector<std::size_t>> element_node_tags;

        gmsh::model::mesh::getElements(element_types, element_tags, element_node_tags);

        // Find triangular elements (type 2)
        std::vector<std::size_t> triangular_elements;
        std::vector<std::size_t> triangular_nodes;

        for (size_t i = 0; i < element_types.size(); ++i)
        {
            if (element_types[i] == 2)
            { // Triangular elements
                triangular_elements = element_tags[i];
                triangular_nodes = element_node_tags[i];
                break;
            }
        }

        size_t num_triangles = triangular_elements.size();
        std::cout << "Number of triangular elements: " << num_triangles << std::endl;

        // 3. Read physical groups and get burner nodes
        std::vector<std::pair<int, int>> physical_groups;
        gmsh::model::getPhysicalGroups(physical_groups);

        std::vector<std::size_t> burner1_nodes, burner2_nodes, burner3_nodes, burner4_nodes;
        std::vector<std::size_t> stove_body_nodes;

        for (const auto &pg : physical_groups)
        {
            int dim = pg.first;
            int tag = pg.second;

            if (dim == 2)
            { // 2D physical groups
                std::string name;
                gmsh::model::getPhysicalName(dim, tag, name);

                std::vector<int> entities;
                gmsh::model::getEntitiesForPhysicalGroup(dim, tag, entities);

                for (int entity : entities)
                {
                    std::vector<std::size_t> nodes;
                    gmsh::model::mesh::getNodes(nodes, node_coords, parametric_coords, dim, entity);

                    if (tag == 100)
                    { // Stove Body
                        stove_body_nodes.insert(stove_body_nodes.end(), nodes.begin(), nodes.end());
                    }
                    else if (tag == 101)
                    { // Burner 1
                        burner1_nodes.insert(burner1_nodes.end(), nodes.begin(), nodes.end());
                    }
                    else if (tag == 102)
                    { // Burner 2
                        burner2_nodes.insert(burner2_nodes.end(), nodes.begin(), nodes.end());
                    }
                    else if (tag == 103)
                    { // Burner 3
                        burner3_nodes.insert(burner3_nodes.end(), nodes.begin(), nodes.end());
                    }
                    else if (tag == 104)
                    { // Burner 4
                        burner4_nodes.insert(burner4_nodes.end(), nodes.begin(), nodes.end());
                    }
                }
            }
        }

        // Remove duplicates from burner node lists
        auto remove_duplicates = [](std::vector<std::size_t> &vec)
        {
            std::sort(vec.begin(), vec.end());
            vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
        };

        remove_duplicates(stove_body_nodes);
        remove_duplicates(burner1_nodes);
        remove_duplicates(burner2_nodes);
        remove_duplicates(burner3_nodes);
        remove_duplicates(burner4_nodes);

        std::cout << "Stove Body nodes: " << stove_body_nodes.size() << std::endl;
        std::cout << "Burner 1 nodes: " << burner1_nodes.size() << std::endl;
        std::cout << "Burner 2 nodes: " << burner2_nodes.size() << std::endl;
        std::cout << "Burner 3 nodes: " << burner3_nodes.size() << std::endl;
        std::cout << "Burner 4 nodes: " << burner4_nodes.size() << std::endl;

        // 4. FEM Formulation and Matrix Assembly
        std::cout << "Assembling global matrices..." << std::endl;

        // Create sparse matrices
        Eigen::SparseMatrix<double> K(num_nodes, num_nodes);
        Eigen::SparseMatrix<double> M(num_nodes, num_nodes);

        // Use triplet list for efficient sparse matrix construction
        std::vector<Eigen::Triplet<double>> K_triplets, M_triplets;

        // Loop through triangular elements
        for (size_t elem = 0; elem < num_triangles; ++elem)
        {
            // Get element nodes (3 nodes per triangle)
            std::size_t n1 = triangular_nodes[3 * elem];
            std::size_t n2 = triangular_nodes[3 * elem + 1];
            std::size_t n3 = triangular_nodes[3 * elem + 2];

            // Get node indices
            size_t i1 = node_tag_to_index[n1];
            size_t i2 = node_tag_to_index[n2];
            size_t i3 = node_tag_to_index[n3];

            // Get node coordinates
            double coords[6] = {
                node_coords[3 * i1], node_coords[3 * i1 + 1], // x1, y1
                node_coords[3 * i2], node_coords[3 * i2 + 1], // x2, y2
                node_coords[3 * i3], node_coords[3 * i3 + 1]  // x3, y3
            };

            // Determine material property based on element location
            double alpha_elem = is_burner_element(coords) ? alpha_burner : alpha_stove;

            // Compute local matrices
            Eigen::Matrix3d k_local = compute_local_stiffness(coords, alpha_elem);
            Eigen::Matrix3d m_local = compute_local_mass(coords);

            // Assembly into global matrices
            std::vector<size_t> indices = {i1, i2, i3};

            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    K_triplets.emplace_back(indices[i], indices[j], k_local(i, j));
                    M_triplets.emplace_back(indices[i], indices[j], m_local(i, j));
                }
            }
        }

        // Build sparse matrices from triplets
        K.setFromTriplets(K_triplets.begin(), K_triplets.end());
        M.setFromTriplets(M_triplets.begin(), M_triplets.end());

        std::cout << "Global matrices assembled successfully" << std::endl;

        // 5. Time-stepping and solving
        std::cout << "Starting time-stepping simulation..." << std::endl;

        // Compute system matrix A = (M/dt) + K
        Eigen::SparseMatrix<double> A = (M / dt) + K;

        // Initialize temperature vector
        Eigen::VectorXd T_old = Eigen::VectorXd::Constant(num_nodes, initial_temp);
        Eigen::VectorXd T_new(num_nodes);

        // Setup solver
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);

        if (solver.info() != Eigen::Success)
        {
            std::cerr << "Error: Matrix decomposition failed!" << std::endl;
            gmsh::finalize();
            return -1;
        }

        // Time-stepping loop
        int step = 0;
        int output_count = 0;

        for (double t = 0; t <= total_time; t += dt, ++step)
        {
            // Construct right-hand side: b = (M/dt) * T_old
            Eigen::VectorXd b = (M / dt) * T_old;

            // Create heat source vector
            Eigen::VectorXd F = Eigen::VectorXd::Zero(num_nodes);

            // Apply heat source to center points of all 4 burners
            std::vector<std::pair<double, double>> burner_centers = {
                {3.0, 1.0}, // Burner 1
                {7.0, 1.0}, // Burner 2
                {3.0, 5.0}, // Burner 3
                {7.0, 5.0}  // Burner 4
            };

            // Calculate realistic heat source (scaled properly for finite elements)
            double heat_power_per_burner = 2000.0; // 2kW per burner (realistic for stove)
            double burner_radius = 0.25;
            double burner_area = M_PI * burner_radius * burner_radius;

            for (const auto &center : burner_centers)
            {
                // Count nodes in burner region first
                int nodes_in_burner = 0;
                for (size_t i = 0; i < num_nodes; ++i)
                {
                    double dx = node_coords[3 * i] - center.first;
                    double dy = node_coords[3 * i + 1] - center.second;
                    double distance = std::sqrt(dx * dx + dy * dy);

                    if (distance <= burner_radius)
                    {
                        nodes_in_burner++;
                    }
                }

                if (nodes_in_burner > 0)
                {
                    double heat_per_node = heat_power_per_burner / nodes_in_burner; // Distribute power among nodes

                    // Apply heat source in a circular region around each burner center
                    for (size_t i = 0; i < num_nodes; ++i)
                    {
                        double dx = node_coords[3 * i] - center.first;
                        double dy = node_coords[3 * i + 1] - center.second;
                        double distance = std::sqrt(dx * dx + dy * dy);

                        if (distance <= burner_radius)
                        {
                            // Apply heat source with Gaussian distribution
                            double heat_factor = std::exp(-distance * distance / (2.0 * 0.1 * 0.1));

                            // Scale down significantly to get realistic temperatures
                            F(i) += heat_per_node * heat_factor * 0.0001; // Much smaller scaling factor
                        }
                    }
                }
            }

            // Add source to right-hand side
            b += F;

            // Solve linear system A * T_new = b
            T_new = solver.solve(b);

            if (solver.info() != Eigen::Success)
            {
                std::cerr << "Error: Linear system solve failed at time " << t << std::endl;
                break;
            }

            // Update temperature
            T_old = T_new;

            // Output at specified intervals
            if (step % output_interval == 0)
            {
                std::string filename = "temperature_" + std::to_string(output_count) + ".pos";
                std::vector<double> temps(T_new.data(), T_new.data() + T_new.size());
                write_pos_file(filename, node_coords, temps);
                output_count++;

                std::cout << "Time: " << std::fixed << std::setprecision(1) << t
                          << ", Max temp: " << std::setprecision(2) << T_new.maxCoeff()
                          << ", Min temp: " << T_new.minCoeff() << std::endl;
            }
        }

        std::cout << "Simulation completed successfully!" << std::endl;
        std::cout << "Generated " << output_count << " output files" << std::endl;

        // Finalize Gmsh
        gmsh::finalize();
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        gmsh::finalize();
        return -1;
    }

    return 0;
}