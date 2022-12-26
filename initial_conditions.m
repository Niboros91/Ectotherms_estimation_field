function x = initial_conditions(Field, Year,N_stages)
    x = zeros(N_stages,1); %Array of variables

    if contains(Year,'2017')
        if contains (Field,'Montelibretti')
            %Initial conditions
            x(1) = 10000000; % eggs
            x(7) = 0; % adult males
            x(8) = 10000000; % adult mated females
        end
    end
     if contains(Year,'2018')
        if contains (Field,'Montelibretti')
             %Initial conditions
            x(1) = 25; % eggs
            x(7) = 0; % adult males
            x(8) = 14; % adult mated females    
        else
              %Initial conditions
            x(1) = 10; % eggs
            x(7) = 0; % adult males
            x(8) = 10; % adult mated females
    
        end
     end
     if contains(Year,'2019')
        if contains (Field,'Montelibretti')
            %Initial conditions
            x(1) = 0; % eggs
            x(7) = 0; % adult males
            x(8) = 2650000; % adult mated females
        else
            %Initial conditions
            x(1) = 0; % eggs
            x(7) = 0; % adult males
            x(8) = 2650000; % adult mated females
        end
    end

end