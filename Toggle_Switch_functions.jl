#Julia functions used to solve numerically the stochastic dynamics system

function mono_togswich_array(_V_rand,_A,_KX,_KY,_gamma,_dt)
     #toogle switch WITHOUT noise
    dX = ( (_A./(_KX.+_V_rand[:,2])).* (_KX.*_V_rand[:,1])./(_KX.+_V_rand[:,1]) .- _V_rand[:,1] .*_gamma ) .*_dt;
    dY = ( (_A./(_KY.+_V_rand[:,1])).* (_KY.*_V_rand[:,2])./(_KY.+_V_rand[:,2]) .- _V_rand[:,2] .*_gamma ) .*_dt;


     d_V=[dX dY];

    return d_V
end

function mono_togswich_stoch_array(_V_rand,_A,_KX,_KY,_gamma,_dt)
     #toogle switch WITH noise
    dX = ( (_A./(_KY.+_V_rand[:,2])).* (_KY.*_V_rand[:,1])./(_KX.+_V_rand[:,1])  .- _V_rand[:,1] .*_gamma ) .*_dt.+(mu_x .*sqrt(dt) .*rand(Normal(),size(_V_rand)[1]));
    dY = ( (_A./(_KX.+_V_rand[:,1])).* (_KX.*_V_rand[:,2])./(_KY.+_V_rand[:,2]) .- _V_rand[:,2] .*_gamma ) .*_dt.+(mu_y .*sqrt(dt) .*rand(Normal(),size(_V_rand)[1]));


     d_V=[dX dY];

    return d_V
end


function mono_togswich_stoch_array_asym(_V_rand,_A,_KX,_KY,_gamma,_dt,_delta_x,_delta_y)

    dX = ( (_A./(_KY.+_V_rand[:,2])).* (_KY.*_V_rand[:,1])./(_KX.+_V_rand[:,1])  .- _V_rand[:,1] .*_gamma .+_delta_x) .*_dt.+(mu_x .*sqrt(dt) .*rand(Normal(),size(_V_rand)[1]));
    dY = ( (_A./(_KX.+_V_rand[:,1])).* (_KX.*_V_rand[:,2])./(_KY.+_V_rand[:,2]) .- _V_rand[:,2] .*_gamma .+_delta_y) .*_dt.+(mu_y .*sqrt(dt) .*rand(Normal(),size(_V_rand)[1]));


     d_V=[dX dY];

    return d_V
end


function mono_togswich_stoch_array_abs(_V_rand,_A,_KX,_KY,_gamma,_dt)

    dX = ( (_A./(_KY.+abs.(_V_rand[:,2]))).* (_KY.*abs.(_V_rand[:,1]))./(_KX.+abs.(_V_rand[:,1]))  .- abs.(_V_rand[:,1]) .*_gamma ) .*_dt.+(mu_x .*sqrt(dt) .*rand(Normal(),size(_V_rand)[1]));
    dY = ( (_A./(_KX.+abs.(_V_rand[:,1]))).* (_KX.*abs.(_V_rand[:,2]))./(_KY.+abs.(_V_rand[:,2])) .- abs.(_V_rand[:,2]) .*_gamma ) .*_dt.+(mu_y .*sqrt(dt) .*rand(Normal(),size(_V_rand)[1]));


     d_V=[dX dY];

    return d_V
end

function outputMonoTS_matrix!( _ns::Float64, _is::IOStream, _t::Float64, _preT::Float64, _ds::Float64)
    #function to print results at each time step simulation 
    if _t > _preT
    tat = _t - _preT
        if (_ns * _ds)  < tat
            #print(_is,tat , "\t")
            print(_is,_t , "\t")
            #for air = 1 : size(X_rand)[1] # loop over array rows
            for air = 1 : size(Mono_TS_matrix)[1] #new line over matrix

                    print(_is, Mono_TS_matrix[air,1] , "\t",Mono_TS_matrix[air,2], "\t")

            end
            print(_is, "\n") #end of array output
            _ns += 1
        end #do
    end # end write out loop
    return _ns
end

#here the function for the concentration in the new coordinate system (x^,y^)

function outputMonoTS_XY_bar_matrix!( _ns::Float64, _is::IOStream, _t::Float64, _preT::Float64, _ds::Float64)
    if _t > _preT
    tat = _t - _preT
        if (_ns * _ds)  < tat
            #print(_is,tat , "\t")
            print(_is,_t , "\t")
            #for air = 1 : size(X_rand)[1] # loop over array rows
            for air = 1 : size(Mono_TS_matrix)[1] #new line over matrix
                X_new=sin(theta) .*(Mono_TS_matrix[air,1]-X0) .-cos(theta).*(Mono_TS_matrix[air,2]-Y0);
                Y_new=cos(theta) .*(Mono_TS_matrix[air,1]-X0) .+sin(theta).*(Mono_TS_matrix[air,2]-Y0);
                print(_is,X_new , "\t",Y_new, "\t")

            end
            print(_is, "\n") #end of array output
            _ns += 1
        end #do
    end # end write out loop
    return _ns
end


#defintion of the function that writes the variance at each time step
function outputVariance_matrix!( _ns::Float64, _is::IOStream, _t::Float64, _preT::Float64, _ds::Float64)
    if _t > _preT
    tat = _t - _preT
        if (_ns * _ds)  < tat
            print(_is,tat , "\t")
            X_ensemble=var(Mono_TS_matrix[:,1])
            Y_ensemble=var(Mono_TS_matrix[:,2])
             print(_is, "\t",X_ensemble , "\t",Y_ensemble, "\t")
             print(_is, "\n") #end of array output for each time step go at the head of the next line
            _ns += 1
        end
    end #do
   return _ns
end


function outputVariance_Mean_matrix!( _ns::Float64, _is::IOStream, _t::Float64, _preT::Float64, _ds::Float64)

    if _t > _preT
    tat = _t - _preT
        if (_ns * _ds)  < tat
            print(_is,tat , "\t")
            X_ensemble=var(Mono_TS_matrix[:,1]; mean=X0)
            Y_ensemble=var(Mono_TS_matrix[:,2]; mean=Y0)
             print(_is, "\t",X_ensemble , "\t",Y_ensemble, "\t")
             print(_is, "\n") #end of array output for each time step go at the head of the next line
            _ns += 1
        end
    end #do
   return _ns
end


function outputVariance_Mean_Norm_matrix!( _ns::Float64, _is::IOStream, _t::Float64, _preT::Float64, _ds::Float64)

    if _t > _preT
    tat = _t - _preT
        if (_ns * _ds)  < tat
            print(_is,tat , "\t")
            X_ensemble=var(Mono_TS_matrix[:,1]; mean=X0)
            Y_ensemble=var(Mono_TS_matrix[:,2]; mean=Y0)
            Variance_Norm=X_ensemble+Y_ensemble
             print(_is, Variance_Norm, "\t")
             print(_is, "\n") #end of array output for each time step go at the head of the next line
            _ns += 1
        end
    end #do
   return _ns
end

#this function evaluates the CV_X and CV_Y
function outputCV_matrix!( _ns::Float64, _is::IOStream, _t::Float64, _preT::Float64, _ds::Float64,_X0::Float64,_Y0::Float64)
    if _t > _preT
    tat = _t - _preT
        if (_ns * _ds)  < tat
            print(_is,tat , "\t")
            sigma_x=std(Mono_TS_matrix[:,1])
            mu_x=abs(mean(Mono_TS_matrix[:,1])-_X0)
            sigma_y=std(Mono_TS_matrix[:,2])
            mu_y=abs(mean(Mono_TS_matrix[:,2])-_Y0)
             print(_is, "\t",sigma_x,"\t",mu_x,"\t",sigma_y,"\t",mu_y,"\t")
             print(_is, "\n") #end of array output for each time step go at the head of the next line
            _ns += 1
        end
    end #do
   return _ns
end

#this function evaluates the CV_X' (induction line) and CV_Y'
function outputCV_inductionline_matrix!( _ns::Float64, _is::IOStream, _t::Float64, _preT::Float64, _ds::Float64,_alpha::Float64)
    if _t > _preT
    tat = _t - _preT
        if (_ns * _ds)  < tat
            print(_is,tat , "\t")
            X_new=cos(_alpha) .* Mono_TS_matrix[:,1].- sin(_alpha).* Mono_TS_matrix[:,2] .+ Kmax .*sin(_alpha) .-KY .*sin(_alpha);
            Y_new=sin(_alpha) .* Mono_TS_matrix[:,1].+ cos(_alpha).* Mono_TS_matrix[:,2] .- Kmax .*cos(_alpha) .+KY .*cos(_alpha);
            X_CV=std(X_new)/mean(X_new)
            Y_CV=std(Y_new)/mean(Y_new)
             print(_is, "\t",X_CV , "\t",Y_CV, "\t")
             print(_is, "\n") #end of array output for each time step go at the head of the next line
            _ns += 1
        end
    end #do
   return _ns
end


function outputVariance_Mean_ontheline_matrix!( _ns::Float64, _is::IOStream, _t::Float64, _preT::Float64, _ds::Float64,_alpha::Float64)

    if _t > _preT
    tat = _t - _preT
        if (_ns * _ds)  < tat
            print(_is,tat , "\t")
            X_new=cos(_alpha) .* Mono_TS_matrix[:,1].- sin(_alpha).* Mono_TS_matrix[:,2] .+ Kmax .*sin(_alpha) .-KY .*sin(_alpha);
            Y_new=sin(_alpha) .* Mono_TS_matrix[:,1].+ cos(_alpha).* Mono_TS_matrix[:,2] .- Kmax .*cos(_alpha) .+KY .*cos(_alpha);
            X_ensemble=var(X_new; mean=0)
            Y_ensemble=var(Y_new; mean=0)
             print(_is, "\t",X_ensemble , "\t",Y_ensemble, "\t")
             print(_is, "\n") #end of array output for each time step go at the head of the next line
            _ns += 1
        end
    end #do
   return _ns
end

function outputMean_ontheline_matrix!( _ns::Float64, _is::IOStream, _t::Float64, _preT::Float64, _ds::Float64,_alpha::Float64)

    if _t > _preT
    tat = _t - _preT
        if (_ns * _ds)  < tat
            print(_is,tat , "\t")
            X_new=cos(_alpha) .* Mono_TS_matrix[:,1].- sin(_alpha).* Mono_TS_matrix[:,2] .+ Kmax .*sin(_alpha) .-KY .*sin(_alpha);
            Y_new=sin(_alpha) .* Mono_TS_matrix[:,1].+ cos(_alpha).* Mono_TS_matrix[:,2] .- Kmax .*cos(_alpha) .+KY .*cos(_alpha);
            X_ensemble=mean(X_new)
            Y_ensemble=mean(Y_new)
             print(_is, "\t",X_ensemble , "\t",Y_ensemble, "\t")
             print(_is, "\n") #end of array output for each time step go at the head of the next line
            _ns += 1
        end
    end #do
   return _ns
end




function outputVariance_Mean_ontheSM_matrix!( _ns::Float64, _is::IOStream, _t::Float64, _preT::Float64, _ds::Float64,_theta::Float64)

    if _t > _preT
    tat = _t - _preT
        if (_ns * _ds)  < tat
            print(_is,tat , "\t")
            X_new=sin(_theta).*(Mono_TS_matrix[:,1] .-X0) .- cos(_theta).*(Mono_TS_matrix[:,2] .- Y0);
            Y_new= cos(_theta) .*(Mono_TS_matrix[:,1] .-X0) .+sin(_theta).*(Mono_TS_matrix[:,2] .- Y0);
            X_ensemble=var(X_new; mean=0)
            Y_ensemble=var(Y_new; mean=0)
             print(_is, "\t",X_ensemble , "\t",Y_ensemble, "\t")
             print(_is, "\n") #end of array output for each time step go at the head of the next line
            _ns += 1
        end
    end #do
   return _ns
end
