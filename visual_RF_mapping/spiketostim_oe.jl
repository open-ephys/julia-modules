#function spiketostim()

fprev = zeros(5001,1);
xprev = zeros(2,5001);



global fprev;
global xprev;
global buffer_count 	=0;	
global nSpikes		=0;
global iter		=0;
global flag_accum = 1;

include("stimBayes.jl");
include("makeSobol.jl");
include("OptProbBayes.jl");
include("EIDMat.jl");
include("kernelD.jl");

#for iter = 1:22
#	nSpikes = rand(100,1);
#	Im = stimBayes(nSpikes, iter);
#end


function oe_process!(data,image)

	global fprev;
	global xprev;
	global buffer_count;	
	global nSpikes;
	global iter;
	global flag_accum::Atomic{Int};

	if (flag_accum==0)
		# do nothing
	else
		buffer_count = buffer_count+1;
	
		if (buffer_count >2 && buffer_count<=30)			
			if (buffer_count > 10)
				image[1:30*30] = 0.5 * ones(1,900);
			else
				nSpikes = nSpikes+sum(data);
			end		
		elseif (buffer_count>30)
			iter	   = iter+1;
			flag_accum = 0;

			# compute next image from previous responses
			image[1:30*30] = stimBayes(nSpikes, iter);		
		
			buffer_count = 0;
	
			data[1] = nSpikes;
			data[2] = xprev[1,iter];
			data[3] = xprev[2,iter];

			println("iter: ",iter);
			nSpikes=0;
			#@printf "nSpikes: %d bufferCount: %f image %f" nSpikes buffer_count image[1];

			flag_accum   = 1;

		end
	end	


	

	return data, image;
end;

