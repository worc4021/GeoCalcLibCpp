classdef GeoCalcLibTest < matlab.unittest.TestCase
    
    properties (Access = private)
        didAddPath (1,1) logical = false
    end

    methods (TestClassSetup)
        function addPath(testCase)
            if isempty(which('vertexEnumeration'))
                testCase.didAddPath = true;
                addpath(fullfile(fileparts(mfilename('fullpath')),'build'));
            end
        end
    end

    methods (TestClassTeardown)
        function removePath(testCase)
            if testCase.didAddPath
                rmpath(fullfile(fileparts(mfilename('fullpath')),'build'));
            end
        end
    end

    methods(Test)
        function testVertexEnumeration(testCase)
            
            n = 3;
            A = [eye(n);-eye(floor(n/2),n)];
            b = ones(n+floor(n/2),1);

            [V,R] = vertexEnumeration(A,b);

            Vref = ones(n-1,n);
            Vref(1,1) = -1;

            Rref = [zeros(n-1,1),-eye(n-1)];

            for v = Vref'
                testCase.verifyTrue(any(all(V==v',2)),"Vertices do not match expected ones.");
            end

            for r = Rref'
                testCase.verifyTrue(any(all(R==r',2)),"Rays do not match expected ones.");
            end
        end


        function testFacetEnumeration(testCase)
            
            V = [-1,0;...
                0,1;...
                0,-1];
            R = [1,0];

            [A,b] = facetEnumeration(V,R);

            Aref = [-1,1;...
                    -1,-1;...
                    0,-1;...
                    0,1];
            bref = ones(4,1);

            tmp = [bref,Aref];
            for hyperplane = tmp'
                testCase.verifyTrue(any(all([b,A]==hyperplane',2)),"Hyperplanes do not match expected ones.");
            end
        end
    end
    
end