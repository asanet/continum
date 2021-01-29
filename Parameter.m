classdef Parameter < handle
        
    properties (Access = public)
        tag
        value = 0
        lower = -realmax
        upper = realmax
        fixed = false
    end
    
    properties (Access = private)
        type = 'par'
    end
       
    methods (Access = public)
        function obj = Parameter(tag, val, lb, ub, fixed)
            if nargin < 1
                error('Parameter needs a "tag".')
            end
            
            obj.tag = tag;
            
            if nargin >= 2 && ~isempty(val)
                obj.value = val;
            end
            if nargin >= 3 && ~isempty(lb)
                obj.lower = lb;
            end
            if nargin >= 4 && ~isempty(ub)
                obj.upper = ub;
            end
            if nargin == 5 && ~isempty(fixed)
                obj.fixed = fixed;
            end                
        end
                
        function Show(obj)

            if obj.fixed
                fx = 'yes';
            else
                fx = 'no';
            end
            fprintf('\n\t Class: %s\n', mfilename)
            fprintf('\t Tag:   %s\n', obj.tag)
            fprintf('\t Value: %2.4e\n', obj.value)
            fprintf('\t Lower: %2.4e\n', obj.lower)
            fprintf('\t Upper: %2.4e\n', obj.upper)
            fprintf('\t Fixed: %s\n', fx)
        end
        
        function type = getType(obj)
            type = obj.type;
        end
    end
end