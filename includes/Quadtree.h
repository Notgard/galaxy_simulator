#pragma once

#include <iostream>
#include <array>
#include <memory>
#include <vector>
#include <algorithm>
#include <stack>
#include <cmath>
#include <cassert>
#include <functional>

#include <iomanip>
#include <limits>
#include <numbers>

#include "Config.h"
#include "Box.h"

#include "Particle.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

// Quadtree implementation :
//  stores the values in a tree structure, where each node has a maximum of 4 children
//  each node has a maximum of QUADTREE_MAX_VALUES values
//  if a node has more than QUADTREE_MAX_VALUES values, it is split into 4 children
//  each child has a box that is a quarter of the parent box
//  the values are stored in the leaf nodes
namespace quadtree
{
    template <typename T, typename GetBox>
    class Quadtree
    {
        static_assert(std::is_convertible_v<std::invoke_result_t<GetBox, const T &>, Box<float>>,
                      "GetBox must be a callable of signature Box<float>(const T&)");

    public:
        struct Node
        {
            // Node(const Box<float> box) : box(box) {}

            // Box<float> box;
            std::array<std::unique_ptr<Node>, QUADTREE_MAX_VALUES> children;
            std::vector<T> values = std::vector<T>();
            Vector2<double> centerOfMass = {0.0, 0.0};
            double totalMass = 0.0;
        };

        Quadtree() = default;

        Quadtree(const Box<float> &box, const GetBox &getBox = GetBox()) : mBox(box), mRoot(std::make_unique<Node>()), mGetBox(getBox) {}

        void add(T value)
        {
            this->mValues.push_back(value);
            add(mRoot.get(), 0, mBox, value);
        }

        void remove(T value) { remove(mRoot.get(), mBox, value); }

        std::vector<T> query(const Box<float> &box) const
        {
            std::vector<T> values;
            query(mRoot.get(), mBox, box, values);
            return values;
        }

        void update_barnes_hut_forces(double dtime)
        {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for (size_t i = 0; i < mValues.size(); i++)
            {
                T value = mValues[i];
                Particle *particle = static_cast<Particle *>(value);
                if (particle->is_sun)
                    continue;
                compute_forces(mRoot.get(), mBox, value, dtime);
            }
        }

        void update_barnes_hut_subtree_forces(Node *node, const Box<float> &box, T value, double dtime)
        {
            Particle *particle = static_cast<Particle *>(value);
            if (particle->is_sun)
                return;
            compute_forces(node, box, value, dtime);
        }

        void update_individual_barnes_hut_forces(T value, double dtime)
        {
            Particle *particle = static_cast<Particle *>(value);
            if (particle->is_sun)
                return;
            compute_forces(mRoot.get(), mBox, value, dtime);
        }

        std::array<Node *, QUADTREE_MAX_VALUES> getChildren(Node *node)
        {
            std::array<Node *, QUADTREE_MAX_VALUES> children;
            for (int i = 0; i < QUADTREE_MAX_VALUES; ++i)
            {
                children[i] = node->children[i].get();
            }
            return children;
        }

        void update_tree_masses()
        {
            update_masses(mRoot.get(), mBox);
        }

        Box<float> getBox() const { return mBox; }

        bool isLeaf(const Node *node) const { return !node->children[0]; }

        const Node *getRoot() const { return mRoot.get(); }

        void printTree() const
        {
            printNode(mRoot.get(), mBox, 0);
        }

        Box<float> computeBox(const Box<float> &box, int i) const
        {
            auto origin = box.getTopLeft();
            auto size = box.getSize() * 0.5f;
            switch (i)
            {
            case NORTH_EAST:
                return {{origin.x + size.x, origin.y}, size};
            case SOUTH_EAST:
                return {{origin.x + size.x, origin.y + size.y}, size};
            case SOUTH_WEST:
                return {{origin.x, origin.y + size.y}, size};
            case NORTH_WEST:
                return {origin, size};
            default:
                return {};
            }
        }

    private:
        static constexpr size_t q_threshold = QUADTREE_MAX_VALUES;
        static constexpr size_t q_max_depth = QUADTREE_MAX_DEPTH;

        Box<float> mBox;
        std::unique_ptr<Node> mRoot;
        GetBox mGetBox;
        std::vector<T> mValues;

        int getQuadrant(const Box<float> &nodeBox, const Box<float> &valueBox) const
        {
            auto center = nodeBox.getCenter();
            bool left = valueBox.left < center.x && valueBox.getRight() < center.x;
            bool right = valueBox.left >= center.x;
            bool top = valueBox.top < center.y && valueBox.getBottom() < center.y;
            bool bottom = valueBox.top >= center.y;

            if (left && top)
                return NORTH_WEST;
            if (right && top)
                return NORTH_EAST;
            if (left && bottom)
                return SOUTH_WEST;
            if (right && bottom)
                return SOUTH_EAST;

            return -1; // Particle doesn't fit entirely into one quadrant
        }

        void add(Node *node, size_t depth, const Box<float> &box, T value)
        {
            assert(node != nullptr);
            Box<float> valueBox = mGetBox(value);
            assert(box.contains(valueBox));

            // If the node is a leaf and already contains a particle, we must split
            if (isLeaf(node) && !node->values.empty())
            {
                split(node, box); // Split immediately if there is already a particle

                // Move existing particles into correct child nodes
                std::vector<T> oldValues = std::move(node->values);
                node->values.clear(); // Clear parent node after splitting

                for (T oldValue : oldValues)
                {
                    int i = getQuadrant(box, mGetBox(oldValue));
                    if (i != -1)
                    {
                        Box<float> childBox = computeBox(box, i);
                        add(node->children[i].get(), depth + 1, childBox, oldValue);
                    }
                }
            }

            // Now, insert the new value into the correct child node
            if (!isLeaf(node))
            {
                int i = getQuadrant(box, mGetBox(value));
                if (i != -1)
                {
                    Box<float> childBox = computeBox(box, i);
                    add(node->children[i].get(), depth + 1, childBox, value);
                }
                else
                {
                    node->values.push_back(value);
                }
            }
            else
            {
                node->values.push_back(value); // This should only happen if the node is truly empty
            }

            // update_masses(node, box);
        }

        void split(Node *node, const Box<float> &box)
        {
            for (int i = 0; i < QUADTREE_MAX_VALUES; ++i)
            {
                node->children[i] = std::make_unique<Node>();
            }

            // Move existing values into child nodes
            std::vector<T> newValues;
            for (T value : node->values)
            {
                int i = getQuadrant(box, mGetBox(value));
                if (i != -1)
                {
                    Box<float> childBox = computeBox(box, i);
                    add(node->children[i].get(), 0, childBox, value); // Insert properly into the child
                }
                else
                {
                    newValues.push_back(value);
                }
            }
            node->values = std::move(newValues);
        }

        bool remove(Node *node, const Box<float> &box, T value)
        {
            assert(node != nullptr);
            assert(box.contains(mGetBox(value)));
            if (isLeaf(node))
            {
                auto it = std::find(node->values.begin(), node->values.end(), value);
                if (it != node->values.end())
                    node->values.erase(it);
                return true;
            }
            int i = getQuadrant(box, mGetBox(value));
            if (i != -1)
                return remove(node->children[i].get(), computeBox(box, i), value);
            return false;
        }

        void query(const Node *node, const Box<float> &nodeBox, const Box<float> &queryBox, std::vector<T> &results) const
        {
            for (T value : node->values)
            {
                results.push_back(value);
            }
            if (!isLeaf(node))
            {
                for (int i = 0; i < QUADTREE_MAX_VALUES; ++i)
                {
                    if (queryBox.intersects(computeBox(nodeBox, i)))
                    {
                        query(node->children[i].get(), computeBox(nodeBox, i), queryBox, results);
                    }
                }
            }
        }

        void printNode(const Node *node, const Box<float> &box, int depth) const
        {
            if (!node)
                return;

            // Indentation for better visualization
            std::cout << std::string(depth * 2, ' ') << "Node Depth: " << depth
                      << " | Box: (" << box.left << ", " << box.top << ") "
                      << "-> (" << box.getRight() << ", " << box.getBottom() << ")"
                      << " | Total mass: " << node->totalMass
                      << " | Center of mass: (" << node->centerOfMass.x << ", " << node->centerOfMass.y << ")"
                      << " | is leaf: " << isLeaf(node)
                      << std::endl;

            // Print values in the node
            // std::cout << std::string(depth * 2, ' ') << "Vector size: " << node->values.size() << std::endl;
            std::cout << std::string(depth * 2, ' ') << "  Values: ";
            for (const auto &value : node->values)
            {
                Particle *particle = static_cast<Particle *>(value);
                Box<float> valueBox = mGetBox(value);
                int i = getQuadrant(box, valueBox);
                std::cout << "(" << particle->position.x << ", " << particle->position.y << ") [q " << i << "] "; // outputs particle positions
            }
            std::cout << "\n";

            // Recursively print child nodes if they exist
            for (int i = 0; i < QUADTREE_MAX_VALUES; ++i)
            {
                if (node->children[i])
                {
                    Box<float> childBox = computeBox(box, i);
                    printNode(node->children[i].get(), childBox, depth + 1);
                }
            }
        }

        void update_masses(Node *node, const Box<float> &box)
        {
            if (!node)
                return; // Avoid null checks later

            double sum_mass = 0.0;
            Vector2<double> weighted_center(0.0, 0.0);

            // Add mass from stored particles
            for (const auto &value : node->values)
            {
                Particle *particle = static_cast<Particle *>(value);
                sum_mass += particle->mass;
                weighted_center += particle->position * particle->mass;
            }

            // Sum up the mass from child nodes
            for (int i = 0; i < QUADTREE_MAX_VALUES; ++i)
            {
                Node *child = node->children[i].get();
                if (child) // Check if child exists before calling update_masses
                {
                    update_masses(child, computeBox(box, i)); // Potentially optimize computeBox

                    sum_mass += child->totalMass;
                    weighted_center += child->centerOfMass * child->totalMass;
                }
            }

            node->totalMass = sum_mass;
            if (sum_mass > 0.0)
            {
                node->centerOfMass = weighted_center / sum_mass;
            }
        }

        void compute_forces(Node *node, const Box<float> &box, T value, double dtime)
        {
            assert(node != nullptr);
            /*
                if (node->values.size() < SMALL_GALAXY_THRESHOLD)
                {
                }
             */
            // if internal node, compute s/d quotient to determine if node is far enough
            // if so, compute force on body from node as single body
            // if not, recursively apply Barnes-Hut algorithm to children
            if (isLeaf(node))
            {
                //std::cout << "Computing force on leaf node" << std::endl;
                // std::cout << "Computing force on leaf node" << std::endl;
                for (T v : node->values)
                {
                    compute_force(value, v, dtime);
                }
            }
            else
            {
                //std::cout << "Computing force on internal node" << std::endl;
                // std::cout << "Computing force on internal node" << std::endl;
                //  compute s/d quotient
                double s = box.width;
                Box<float> valueBox = mGetBox(value);
                Vector2<float> boxCenter = valueBox.getCenter();
                Vector2<double> center = {boxCenter.x, boxCenter.y};
                double d = Vector2<double>::distance(center, node->centerOfMass);
/* 
                std::cout << "----------------------------------" << std::endl;
                std::cout << "Box center: " << boxCenter.x << ", " << boxCenter.y << std::endl;
                std::cout << "Node center of mass: " << node->centerOfMass.x << ", " << node->centerOfMass.y << std::endl;
                std::cout << "Box width: " << s << std::endl;
                std::cout << "Distance: " << d << std::endl;
                std::cout << "s/d: " << s / d << std::endl;
                std::cout << "Theta: " << THETA << std::endl;
                std::cout << "----------------------------------" << std::endl;
 */
                if (s / d < THETA)
                {
                    // std::cout << "Node " << node << " is far enough" << std::endl;
                    //   compute force on body from node as single body
                    compute_force(value, node, dtime);
                }
                else
                {
                    // recursively apply Barnes-Hut algorithm to children
                    for (int i = 0; i < QUADTREE_MAX_VALUES; ++i)
                    {
                        if (node->children[i]) // Ensure child node exists before processing
                        {
                            compute_forces(node->children[i].get(), computeBox(box, i), value, dtime);
                        }
                    }
                }
            }
        }

        void compute_force(T value, Node *node, double dtime)
        {
            Particle *particle = static_cast<Particle *>(value);
            if (particle->is_sun)
                return;

            Vector2<double> a_g(0.0);
            Vector2<double> r = particle->position - node->centerOfMass;
            double r_mag = Vector2<double>::length(r);

            if (r_mag > 0)
            {
                //std::cout << "Vector r: " << r.x << ", " << r.y << std::endl;
                double acceleration = -G * node->totalMass / (r_mag * r_mag);
                Vector2<double> r_unit_vector = r / r_mag;
                a_g += acceleration * r_unit_vector;
            }
            //std::cout << "Acceleration: " << a_g.x << ", " << a_g.y << std::endl;

            // Update particle acceleration
            particle->update_acceleration(a_g);
        }

        void compute_force(T v1, T v2, double dtime)
        {
            Particle *p1 = static_cast<Particle *>(v1);
            Particle *p2 = static_cast<Particle *>(v2);
            // Avoid self-interaction
            if (p1->id == p2->id)
            {
                return;
            }
            else if (p1->is_sun && p2->is_sun)
            {
                return;
            }

            Vector2<double> a_g(0.0);
            Vector2<double> r = p1->position - p2->position;
            double r_mag = Vector2<double>::length(r);

            if (r_mag > 0)
            {
                //std::cout << "Vector r: " << r.x << ", " << r.y << std::endl;
                double acceleration = -G * p2->mass / (r_mag * r_mag);
                Vector2<double> r_unit_vector = r / r_mag;
                a_g += acceleration * r_unit_vector;
            }
            //std::cout << "Acceleration: " << a_g.x << ", " << a_g.y << std::endl;

            p1->update_acceleration(a_g);
        }
    };
}