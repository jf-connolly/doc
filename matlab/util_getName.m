function name = util_getName(nameTest, typeTest)

name = struct('ln', '', 'hp', '', 'al', '', 'legend', '');

if strcmp(nameTest, 'total')
    switch typeTest

        case 1, name.ln     = 'Inc';
                name.hp     = 'EoMA+d';
                name.al     = 'multi';
                name.legend = sprintf('ADNPSO');
                name.width  = 100;
                
        case 2, name.ln     = 'Inc';
                name.hp     = 'EoMA+t';
                name.al     = 'dnpso';
                name.legend = sprintf('DNSPO');
                name.width  = 1600;
                
        case 3, name.ln     = 'Inc';
                name.hp     = 'EoMA+t';
                name.al     = 'mopso';
                name.legend = sprintf('MOPSO');
                name.width  = 0;

        case 4, name.ln     = 'Inc';
                name.hp     = 'Hdnc';
                name.al     = 'dnpso';
                name.legend = sprintf('GBEST');
                name.width  = 1600;
                
        case 5, name.ln     = 'Bth';
                name.hp     = '';
                name.al     = 'Knn';
                name.legend = sprintf('$k$NN');
                name.width  = 0;
                
        case 6, name.ln     = 'Bth';
                name.hp     = 'EoFAM+t';
                name.al     = 'pso';
                name.legend = sprintf('$\\mathit{EoFAM}_t^\\mathrm{all-B}$');
                name.width  = 1600;
    end

elseif strcmp(nameTest, 'dnpso')
    switch typeTest

        case 1, name.ln     = 'Inc';
                name.hp     = 'EoLB+d';
                name.al     = 'dnpso';
                name.legend = sprintf('');
                name.width  = 15000;

        case 2, name.ln     = 'Inc';
                name.hp     = 'EoFAM+t';
                name.al     = 'dnpso';
                name.legend = sprintf('');
                name.width  = 15000;
                
        case 3, name.ln     = 'Inc';
                name.hp     = 'Hdnc';
                name.al     = 'dnpso';
                name.legend = sprintf('');
                name.width  = 15000;

        case 4, name.ln     = 'Bth';
                name.hp     = '';
                name.al     = 'Knn';
                name.legend = sprintf('');
                name.width  = 0;
 			
    end
elseif strcmp(nameTest, 'mopso')
    switch typeTest

        case 1, name.ln     = 'Inc';
                name.hp     = 'EoMA+t';
                name.al     = 'mopso';
                name.legend = sprintf('$\\mathit{EoMA}_t^\\mathrm{all}$');

        case 2, name.ln     = 'Inc';
                name.hp     = 'Hdnc';
                name.al     = 'mopso';
                name.legend = sprintf('$\\mathit{EoMA}_t$');
    end
elseif strcmp(nameTest, 'video')
    switch typeTest

        case 1, name.ln     = 'Inc';
                name.hp     = 'EoMA+d';
                name.al     = 'multi';
                name.legend = sprintf('APSO');
                name.width  = 1000;

        case 2, name.ln     = 'Inc';
                name.hp     = 'EoLB+d';
                name.al     = 'dnpso';
                name.legend = sprintf('DNPSO');
                name.width  = 1600;
                
        case 3, name.ln     = 'Inc';
                name.hp     = 'EoMA+t';
                name.al     = 'mopso';
                name.legend = sprintf('MOPSO');
                name.width  = 0;
                
        case 4, name.ln     = 'Bth';
                name.hp     = '';
                name.al     = 'Knn';
                name.legend = sprintf('');
                name.width  = 0;
    end
elseif strcmp(nameTest, 'vid')
    switch typeTest

        case 1, name.ln     = 'Inc';
                name.hp     = 'EoLB+d';
                name.al     = 'dnpso';
                name.legend = sprintf('LBEST$_\\mathrm{+d}$');
                name.width  = 15000;

        case 2, name.ln     = 'Inc';
                name.hp     = 'EoFAM+t';
                name.al     = 'dnpso';
                name.legend = sprintf('SWARM');
                name.width  = 15000;
                
        case 3, name.ln     = 'Inc';
                name.hp     = 'Hdnc';
                name.al     = 'dnpso';
                name.legend = sprintf('GBEST');
                name.width  = 15000;
                
        case 4, name.ln     = '';
                name.hp     = '';
                name.al     = '';
                name.legend = sprintf('PSO$_\\mathrm{B}$');
                name.width  = 0;
                
        case 5, name.ln     = '';
                name.hp     = '';
                name.al     = '';
                name.legend = sprintf('$k$NN');
                name.width  = 0;
    end
else
    fprintf('/*-- Error - no name...\n');
end
