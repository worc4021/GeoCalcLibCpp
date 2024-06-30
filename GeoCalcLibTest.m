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

            [Vall,type,linearities] = vertexEnumeration(A,b);
            V = Vall(type,:);
            R = Vall(~type,:);

            Vref = ones(n-1,n);
            Vref(1,1) = -1;

            Rref = [zeros(n-1,1),-eye(n-1)];

            for v = Vref'
                testCase.verifyTrue(any(all(V==v',2)),"Vertices do not match expected ones.");
            end

            for r = Rref'
                testCase.verifyTrue(any(all(R==r',2)),"Rays do not match expected ones.");
            end

            testCase.verifyEmpty(linearities, "Linearities detected");
        end


        function testFacetEnumeration(testCase)
            
            V = [-1,0;...
                0,1;...
                0,-1];
            R = [1,0];

            data = [V;R];
            isVertex = [true(3,1);false(1,1)];

            [A,b,linearities,vol] = facetEnumeration(data,isVertex);

            allIdx = 1:numel(b);
            Aineq = A(setdiff(allIdx,linearities),:);
            bineq = b(setdiff(allIdx,linearities));

            Aref = [-1,1;...
                    -1,-1;...
                    0,-1;...
                    0,1];
            bref = ones(4,1);

            tmp = [bref,Aref];
            for hyperplane = tmp'
                testCase.verifyTrue(any(all([bineq,Aineq]==hyperplane',2)),"Hyperplanes do not match expected ones.");
            end
            testCase.verifyEmpty(linearities, "Linearities were present");
            
        end

        function testLinearities(testCase)
            V = 3*eye(2);
            [A,b,linearities] = facetEnumeration(V);

            allIdx = 1:numel(b);
            Aineq = A(setdiff(allIdx,linearities),:);
            bineq = b(setdiff(allIdx,linearities));
            Aeq = A(linearities,:);
            beq = b(linearities);

            testCase.verifyEqual(size(Aineq,1),2,"Number of inequalities");
            testCase.verifyLessThanOrEqual(Aineq*(V(1,:)'),bineq,"Inequalities");
            testCase.verifyLessThanOrEqual(Aineq*(V(2,:)'),bineq,"Inequalities");
            testCase.verifyEqual(size(Aeq,1),1, "Number of linearities");
            testCase.verifyEqual(Aeq*(V(1,:)'),beq,"Linearity");
            testCase.verifyEqual(Aeq*(V(2,:)'),beq,"Linearity");
        end

        function testVolume(testCase)
            V = [1,1;
                -1,1;
                1,-1;
                -1,-1];

            
            [~,~,~,refVol] = facetEnumeration(V);

            testCase.verifyEqual(refVol,4,"Volume");

            for i = 1:10
                alpha = rand()*2*pi;
                P = [cos(alpha),sin(alpha);-sin(alpha),cos(alpha)];
                [~,~,~,vol] = facetEnumeration(V*P);
                testCase.verifyLessThanOrEqual(abs(refVol-vol),1e-13,"Volume once rotated");
            end

        end
    end
    
end