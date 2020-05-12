function make(ownGmpClass)
    arguments
        ownGmpClass (1,1) logical = true;
    end

    IEIGEN = '/usr/local/Cellar/eigen/3.3.7/include/eigen3';
    IBUILDINGBLOCKS = '/Users/Manuel/Documents/Development/MexBuildingBlocks';
    ILRS = '/Users/Manuel/Documents/Development/GeoCalcLibCpp/lrslib-070';
    IGMP = '/usr/local/Cellar/gmp/6.2.0';

    obj = ["lrsinterface.cpp","mpq.cpp","mpz.cpp",...
        fullfile(IBUILDINGBLOCKS,'eigenConversion.cpp'),...
        fullfile(ILRS,'lrslib.c'),...
        fullfile(ILRS,'lrsgmp.c')];
    if ownGmpClass
        obj = [obj,"mpqMatrix.cpp"];
    else
        obj = [obj,"gmpMatrix.cpp"];
    end
    
    products = ["facetEnumeration","vertexEnumeration"];

    mexArgs = {'-outdir',pwd(),...,
        ['-I',IEIGEN],...
        ['-I',IBUILDINGBLOCKS],...
        ['-I',ILRS],...
        ['-I',fullfile(IGMP,'include')],...
        '-DGMP',...
        '-DMP',...
        '-lgmp'};
    
    if ismac()
        mexArgs = [mexArgs,...
        'CC=/usr/bin/xcrun  clang',...
        'CXX=/usr/bin/xcrun  clang++',...
        'LD=/usr/bin/xcrun  clang',...
        'LDXX=/usr/bin/xcrun  clang++',...
        'CXXFLAGS=-std=c++17 -Wall'];
    end
    
    for i = 1:numel(obj)
        mex(mexArgs{:},'-c',char(obj(i)))
        [~,object,~] = fileparts(obj(i));
        obj(i) = object + ".o";
    end
    
    
    for i = 1:length(products)
        mex('-output',char(products(i)),...
            mexArgs{:},...
            obj.cellstr{:},...
            [char(products(i)),'.cpp']);
    end
    
    delete('*.o');
    
end