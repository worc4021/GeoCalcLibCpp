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
            utilities::error("At least vertices have to be supplied");
        
        if (!utilities::isnumeric(inputs[0]))
            utilities::error("Vertices have to be numeric.");

        Polytope poly;

        if (inputs.size()>1) {
            if (inputs[0].getDimensions()[1] != inputs[1].getDimensions()[1])
                utilities::error("Vertices and rays have to be in same dimension");
            
            poly.setRays(fromMatrixXd(utilities::eigen::convert(inputs[1])));
        }

        
        poly.setVertices(fromMatrixXd(utilities::eigen::convert(inputs[0])));


        if (outputs.size()>0)
            outputs[0] = utilities::eigen::convert(fromMatrixXq(poly.getA()));

        if (outputs.size()>1)
            outputs[1] = utilities::eigen::convert(fromMatrixXq(poly.getB()));

        if (outputs.size() > 2)
            outputs[2] = utilities::eigen::convert(fromMatrixXq(poly.getAeq()));

        if (outputs.size() > 3)
            outputs[3] = utilities::eigen::convert(fromMatrixXq(poly.getBeq()));

    }
};