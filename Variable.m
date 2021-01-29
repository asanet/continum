classdef Variable < handle
        
    properties (Access = public)
        tag
        value = 0
        lower = -realmax
        upper = realmax
    end
    
    properties (Access = private)
        type = 'var'
    end
       
    methods (Access = public)
        function obj = Variable(tag, val, lb, ub)
            if nargin < 1
                error('Variable needs a "tag".')
            end
            
            obj.tag = tag;
            
            if nargin >= 2 && ~isempty(val)
                obj.value = val;
            end
            if nargin >= 3 && ~isempty(lb)
                obj.lower = lb;
            end
            if nargin == 4 && ~isempty(ub)
                obj.upper = ub;
            end
            
        end
        
        function Show(obj)
            fprintf('\n\t Class: %s\n', mfilename)
            fprintf('\t Tag: %s\n', obj.tag)
            fprintf('\t Value: %2.4e\n', obj.value)
            fprintf('\t Lower: %2.4e\n', obj.lower)
            fprintf('\t Upper: %2.4e\n', obj.upper)
        end
        
        function type = getType(obj)
            type = obj.type;
        end
    end
end