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
        
        if (inputs.size() < 1)
            displayError("At least vertices have to be supplied");
        
        if (!isnumeric(inputs[0]))
            displayError("Vertices have to be numeric.");

        

        Polytope poly;

        if (inputs.size()>1) {
            if (inputs[0].getDimensions()[1] != inputs[1].getDimensions()[1])
                displayError("Vertices and rays have to be in same dimension");
            
            poly.setRays(fromMatrixXd(convert(inputs[1])));
        }

        
        poly.setVertices(fromMatrixXd(convert(inputs[0])));


        if (outputs.size()>0)
            outputs[0] = convert(factory, fromMatrixXq(poly.getA()));

        if (outputs.size()>1)
            outputs[1] = convert(factory, fromMatrixXq(poly.getB()));

    }
};