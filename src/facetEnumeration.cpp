#include "mex.hpp"
#include "mexAdapter.hpp"
#include "eigen/conversions.hpp"
#include "lrsinterface.hpp"

FILE *lrs_ifp;
FILE *lrs_ofp;

class MexFunction : public matlab::mex::Function {

public:
    MexFunction() {
        matlabPtr = getEngine();
    }
    ~MexFunction() = default;
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        
        if (inputs.size() < 1)
            utilities::error("List of vertices has to be supplied");
        
        if (!utilities::isnumeric(inputs[0]))
            utilities::error("Vertices have to be numeric.");

        Polytope poly;

        MatrixXq raysAndVertices = fromMatrixXd(utilities::eigen::convert(inputs[0]));
        std::vector<bool> isVertex(raysAndVertices.rows(), true);

        if (inputs.size()>1) {
            if (inputs[0].getDimensions()[0] != inputs[1].getNumberOfElements())
                utilities::error("Vertices and rays have to be in same dimension");
            matlab::data::TypedArray<bool> isVertexArray = std::move(inputs[1]);
            std::transform(isVertexArray.begin(), isVertexArray.end(), isVertex.begin(),[](auto i) { return i; });
        }

        poly.setVrep(raysAndVertices, isVertex);

        Hrep hrep;
        poly.getHrep(hrep);

        MatrixXq Hh(hrep.rows.size(), hrep.nDim + 1);
        for (std::size_t iRow = 0; iRow < hrep.rows.size(); iRow++) {
            Hh.row(iRow) = hrep.rows[iRow];
        }
        
        Eigen::MatrixXd A = fromMatrixXq(Hh.rightCols(hrep.nDim));
        Eigen::VectorXd b = fromMatrixXq(Hh.col(0));

        if (outputs.size()>0)
            outputs[0] = utilities::eigen::convert(A);

        if (outputs.size()>1)
            outputs[1] = utilities::eigen::convert(b);

        if (outputs.size() > 2){
            matlab::data::ArrayFactory factory;
            matlab::data::TypedArray<double> linearities = factory.createArray<double>({hrep.linearities.size(),1});
            std::transform(hrep.linearities.begin(), hrep.linearities.end(), linearities.begin(), [](auto i) { return static_cast<double>(i + 1); });
            outputs[2] = std::move(linearities);
            if (outputs.size() > 3) {
                outputs[3] = factory.createScalar(poly.volume);
            }
        }


    }
};