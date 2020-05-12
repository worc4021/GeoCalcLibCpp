#include "eigenConversion.hpp"
#include "mexFunctions.hpp"
#include "lrsinterface.hpp"

FILE *lrs_ifp;
FILE *lrs_ofp;

typedef matlab::data::Array mexArray;

class MexFunction : public matlab::mex::Function {
private:
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
    ArrayFactory factory;
    void displayError(std::string errorMessage) {
        matlabPtr->feval(
                        matlab::engine::convertUTF8StringToUTF16String("error"),
                        0, 
                        std::vector<mexArray>({
                            factory.createScalar(errorMessage) 
                            })
                        );
    }

public:
    MexFunction() : matlabPtr(getEngine()) {}
    ~MexFunction() {}
    void operator()(ArgumentList outputs, ArgumentList inputs) {
        
        if (inputs.size() < 2)
            displayError("A,b have to be supplied");
        
        if (!isnumeric(inputs[0]) || !isnumeric(inputs[1]))
            displayError("Data has to be numeric.");

        Polytope poly;

        if (inputs.size() > 2) {
            if (inputs.size() < 4)
                displayError("If Aineq and Aeq are provided four arguments are required");
            
            if (inputs[0].getDimensions()[1] != inputs[2].getDimensions()[1])
                displayError("Aineq and Aeq have to be in the same dimension");
            
            poly.setHrepEq(fromMatrixXd(convert(inputs[2])), fromMatrixXd(convert(inputs[3])));
        }
        

        poly.setHrepIneq(fromMatrixXd(convert(inputs[0])), fromMatrixXd(convert(inputs[1])));

        if (outputs.size()>0)
            outputs[0] = convert(factory, fromMatrixXq(poly.getVertices()));

        if (outputs.size()>1)
            outputs[1] = convert(factory, fromMatrixXq(poly.getRays()));
        
        if (outputs.size()>2)
            outputs[2] = factory.createScalar<double>(poly.getVolume());

    }
};