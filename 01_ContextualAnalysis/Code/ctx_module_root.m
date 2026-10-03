function moduleRoot = ctx_module_root()
% Return the root directory of the contextual-analysis module.
codeDir = string(fileparts(mfilename("fullpath")));
moduleRoot = string(fileparts(codeDir));
end
