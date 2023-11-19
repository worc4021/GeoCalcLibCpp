#include "mex.hpp"
#include "mexAdapter.hpp"
#include "eigen/conversions.hpp"
#include "lrsinterface.hpp"

FILE *lrs_ifp;
FILE *lrs_ofp;

class MexFunction : public matlab::mex::Function
{

public:
    MexFunction()
    {
        matlabPtr = getEngine();
    }
    ~MexFunction() = default;
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs)
    {

        if (inputs.size() < 2)
            utilities::error("A,b have to be supplied");

        if (!utilities::isnumeric(inputs[0]) || !utilities::isnumeric(inputs[1]))
            utilities::error("Data has to be numeric.");

        Polytope poly;

        if (inputs.size() > 2)
        {
            if (inputs.size() < 4)
                utilities::error("If Aineq and Aeq are provided four arguments are required");

            if (inputs[0].getDimensions()[1] != inputs[2].getDimensions()[1])
                utilities::error("Aineq and Aeq have to be in the same dimension");

            poly.setHrepEq(fromMatrixXd(utilities::eigen::convert(inputs[2])), fromMatrixXd(utilities::eigen::convert(inputs[3])));
        }

        poly.setHrepIneq(fromMatrixXd(utilities::eigen::convert(inputs[0])), fromMatrixXd(utilities::eigen::convert(inputs[1])));

        if (outputs.size() > 0)
            outputs[0] = utilities::eigen::convert(fromMatrixXq(poly.getVertices()));

        if (outputs.size() > 1)
            outputs[1] = utilities::eigen::convert(fromMatrixXq(poly.getRays()));

        if (outputs.size() > 2)
        {
            matlab::data::ArrayFactory factory;
            outputs[2] = factory.createScalar<double>(poly.getVolume());
        }
    }
};