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

        std::vector<std::size_t> linearities(0);

        MatrixXq A = fromMatrixXd(utilities::eigen::convert(inputs[0]));
        MatrixXq b = fromMatrixXd(utilities::eigen::convert(inputs[1]));

        if (inputs.size() > 2)
        {
            matlab::data::TypedArray<double> linearitiesArray = std::move(inputs[2]);
            linearities.resize(linearitiesArray.getNumberOfElements());
            std::transform(linearitiesArray.begin(), linearitiesArray.end(), linearities.begin(), [](auto i) { return static_cast<std::size_t>(i) - 1; });
        }
        Vrep vrep;
        poly.setHrep(A, b, linearities);
        poly.getVrep(vrep);

        MatrixXq verticesAndRays(vrep.rows.size(), vrep.nDim);
        
        matlab::data::ArrayFactory factory;
        matlab::data::TypedArray<bool> isVertex = factory.createArray<bool>({vrep.rows.size(),1});
        qType one_t(1L);
        for (std::size_t iRow = 0; iRow < vrep.rows.size(); iRow++)
        {
            verticesAndRays.row(iRow) = vrep.rows[iRow].tail(vrep.nDim);
            isVertex[iRow] = vrep.rows[iRow](0) == one_t;
        }
        
        if (outputs.size() > 0)
            outputs[0] = utilities::eigen::convert(fromMatrixXq(verticesAndRays));

        if (outputs.size() > 1)
            outputs[1] = std::move(isVertex);

        if (outputs.size() > 2) {
            matlab::data::TypedArray<double> linearitiesArray = factory.createArray<double>({vrep.linearities.size(),1});
            std::transform(vrep.linearities.begin(), vrep.linearities.end(), linearitiesArray.begin(), [](auto i) { return static_cast<double>(i+1); });
            outputs[2] = std::move(linearitiesArray);
        }

    }
};