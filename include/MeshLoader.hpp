//
// Created by Enrico on 25/11/2023.
//

#ifndef EIKONEL_TEST_MESHLOADER_H
#define EIKONEL_TEST_MESHLOADER_H

#include <Mesh.h>
#include <VtkParser.hpp>

class MeshLoader {
public:
    static int load(Mesh &mesh, const VtkParser &parser, std::unordered_map<Point,double>& pointsData,
                    std::unordered_map<MeshElement,double>& elementData);

    static int dump(const Mesh &mesh, VtkParser &parser, const std::unordered_map<Point,double> &pointsData,
                    const std::unordered_map<MeshElement,double>& elementData);
};


#endif //EIKONEL_TEST_MESHLOADER_H
