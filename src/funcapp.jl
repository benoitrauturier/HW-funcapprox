module funcapp
	using FastGaussQuadrature
	using ApproxFun
	using PyPlot
	using ApproXD
	using CompEcon
#Defining the mathematical functions used throuout the homework
	function f(x)
		return x .+ 2.0.*x.^2 - exp(-x)
	end
	function f2(x)
		return 1.0 ./ (1 .+ 25 .* x.^2)
	end
	function f3(x)
		return abs(x).^(0.5)
	end


	# use chebyshev to interpolate this:
	function q1(n)
########################INTERPOLATION NODES#####################################
	#Here we set up the interpolation nodes and normalized interpolation nodes as in slide #10
		nodes=zeros(n)
		for i in 1:n
			val=(n-i+0.5)/(n)
			if  (val % 2 != 0.5)&(val % 2 != 1.5)
				nodes[i]=3.0*cos(pi*val)
			else
				nodes[i]=0.0
			end
		end
		#And we compute the value of f at the nodes
		fvalues=f(nodes)
########################INTERPOLATION MATRIX####################################
		Phi=zeros(n,n-1)
		for i in 1:n
			for j in 1:(n-1)
				val=(n-i+0.5)*(j-1.0)/n
				if  (val % 2 != 0.5)&(val % 2 != 1.5)
					Phi[i,j]=cos(pi*val)
				else
					Phi[i,j]=0.0
				end
			end
		end
	#Here we compute c. Since N and J are different we use the least square estimator
		c=(transpose(Phi)*Phi)^(-1)*transpose(Phi)*fvalues
#######################SET UP OF THE INTERPOLATOR###############################
		function interpolator(x)
			#We construct the interpolation matrix for x
			x_11=(2.0*(x+3.0)/6.0)-1.0
			phi_x=zeros(n-1)
			phi_x[1]=1.0
			phi_x[2]=x_11
			for i in 3:(n-1)
				phi_x[i]=2.0*x_11*phi_x[i-1]-phi_x[i-2]
			end
			y=transpose(phi_x)*c
			return y
		end
		return interpolator
	end
		# rest of the question in function runall()
##################################END OF QUESTION 1############################
	function q2(n)
		S=Chebyshev([-3,3])
		x=points(S,n)
		#True value of Chebyshev nodes
		val=f(x)
		#Construction of the interpolator
		interpolator=Fun(ApproxFun.transform(S,val),S)
		return interpolator
	end
	#rest of question in function runall
###################################END OF QUESTION 2############################

	# plot the first 9 basis Chebyshev Polynomial Basisi Fnctions
	function q3()
		#I am pretty sure this is not the fastest way to compute a Chebyshev Polynomial
		#It should be more precise than with the arc cos function (even though the error in that case is of oreder e-16)
		function Chebypol(x,n)
			if n==0
				return 1
			elseif n==1
			 	return x
			else
				return 2.0.*x.*Chebypol(x,n-1)-Chebypol(x,n-2)
			end
		end
		#Ploting the polynomials
		for n in 1:9
			order=n-1
			x=linspace(-1,1,50)
			y=zeros(50)
			for i in 1:50
				y[i]=Chebypol(x[i],order)
			end
			if n==1
				figure("Question 3")
			end
			numplot=330 + n
			subplot(numplot)
			plot(x,y,color="red")
			title("Order $order")
		end
	end
#############################END OF QUESTION 3#################################
	ChebyT(x,deg) = cos(acos(x)*deg)
	unitmap(x,lb,ub) = 2.*(x.-lb)/(ub.-lb) - 1	#[a,b] -> [-1,1]

	type ChebyType
		f::Function # fuction to approximate
		nodes::Union{Vector,LinSpace} # evaluation points
		basis::Matrix # basis evaluated at nodes
		coefs::Vector # estimated coefficients

		deg::Int 	# degree of chebypolynomial
		lb::Float64 # bounds
		ub::Float64

		# constructor
		function ChebyType(_nodes::Union{Vector,LinSpace},_deg,_lb,_ub,_f::Function)
			n = length(_nodes)
			y = _f(_nodes)
			_basis = Float64[ChebyT(unitmap(_nodes[i],_lb,_ub),j) for i=1:n,j=0:_deg]
			_coefs = _basis \ y  # type `?\` to find out more about the backslash operator. depending the args given, it performs a different operation
			# create a ChebyType with those values
			new(_f,_nodes,_basis,_coefs,_deg,_lb,_ub)
		end
	end

	# function to predict points using info stored in ChebyType
	function predict(Ch::ChebyType,x_new)

		true_new = Ch.f(x_new)
		basis_new = Float64[ChebyT(unitmap(x_new[i],Ch.lb,Ch.ub),j) for i=1:length(x_new),j=0:Ch.deg]
		basis_nodes = Float64[ChebyT(unitmap(Ch.nodes[i],Ch.lb,Ch.ub),j) for i=1:length(Ch.nodes),j=0:Ch.deg]
		preds = basis_new * Ch.coefs
		preds_nodes = basis_nodes * Ch.coefs

		return Dict("x"=> x_new,"truth"=>true_new, "preds"=>preds, "preds_nodes" => preds_nodes)
	end

##################################BEGINING OF QUESTION 4a#######################

	function q4a(deg=(5,9,15),lb=-5.0,ub=5.0)
		for i in 1:length(deg)
			#All variables with a 2 are for Chebyshev nodes related result. Without any numbers it is for equidistant interpolation nodes
			#We use linspace for equidistant nodes and we create the cheby nodes
			nodes=linspace(lb,ub,deg[i]+1)
			nodes2=zeros(deg[i]+1)
			for j in 1:(deg[i]+1)
				nodes2[j]=((lb+ub)/2)+((ub-lb)/2)*cos(pi*(deg[i]+1-j+0.5)/(deg[i]+1))
			end
			#we construct the predictor from the class ChebyType
			predictor=ChebyType(nodes,deg[i],lb,ub,f2)
			predictor2=ChebyType(nodes2,deg[i],lb,ub,f2)
			#And compute the results.
			result=predict(predictor, linspace(lb,ub,500))
			result2=predict(predictor2, linspace(lb,ub,500))
			#In the first loop we set up the figure as well as the graph of the true value
			if i==1
				figure("Question 4a")
				subplot(241)
				plot(result["x"],result["truth"], color="red")
				subplot(245)
				plot(result2["x"],result2["truth"], color="red")
			end
			#and at each loop we plot the approximation
			subplot(241 + i)
			plot(result["x"],result["preds"], color="blue")
			title("Uniform")
			subplot(245 + i)
			plot(result2["x"],result2["preds"], color="blue")
			title("Cheby")
		end
		#=
		the pike at x=0 is very poorly represented by the interpolation functions.
		I guess that is the purpose of the exercise.
		The uniform interpolation performs even more poorly
		=#

	end

##################################END OF QUESTION 4a##########################

	function q4b()
		version1=BSpline(13,3,-5,5)
		#The following yield knots more concentrated around zero that to the boundary of the interval
		#=We first use Chebyshev nodes. As they are more concentrated towards the
		boundary we make the negative nodes positive (by adding 5) and the positive
		one negative We keep the nodes at the boundary unchanged=#
		cknots=zeros(13)
		for i in 1:13
			cknots[i]=5*cos(pi*(13-i+0.5)/(13))
			if abs(cknots[i])<0.1 #To avoid to have two points very close to a boundary
				continue
			elseif ((cknots[i]<0.0) & (cknots[i]>-4.9))
				cknots[i]=cknots[i]+5
			elseif ((cknots[i]<4.9) & (cknots[i]>0.0))
				cknots[i]=cknots[i]-5
			elseif cknots[i]<-4.9
				cknots[i]=-5.0
			elseif cknots[i]>4.9
				cknots[i]=5.0
			end
		end
		#We sort the knots in the right order as this is imposed by BSpline function
		cknots=sort(cknots)
		version2=BSpline(cknots,13)
		#defining the equidistant evaluation points
		nevals=linspace(-5,5,65)
		#We compute the base for both group of knots
		basis1 = full(getBasis(collect(nevals),version1))
		basis2 = full(getBasis(collect(nevals),version2))
		#true value
		y=f2(nevals)
		#And computing the coefficients
		coef1 = basis1 \ y
		coef2 = basis2 \ y
		#Computing the estimation for 500 points, for both strategies
		points=linspace(-5,5,500)
		base1=full(getBasis(collect(points),version1))
		base2=full(getBasis(collect(points),version2))
		estimation1= base1 * coef1
		estimation2= base2 * coef2
		truth = f2(points)
		#Ploting the results
		figure("Question 4b")
		subplot(121)
		plot(points,estimation1,color="red")
		plot(points,estimation2,color="green")
		plot(points,truth,color="blue")
		title("true value in blue, equidistant interpolation points in red, concentrated interpolation points in green")
		subplot(122)
		plot(points, truth .- estimation1, color="red")
		plot(points, truth .- estimation2, color="green")
		plot(cknots,zeros(13),"go")
		#=
		Here one sees that concentrating the nodes around the pic reduces the error.
		Not by much in my example, probably because my nodes are not concentrated enough.
		On the countrary it seems to increase the bias at the boundary of the interval.
		=#
	end
###################################END OF QUESTION 4b##########################

	function q5(n=2)
		#Define n as the number of knots =0
		points=collect(linspace(-1,1,500))
		truth=f3(points)
		#We evaluate the points as before
		version1=BSpline(13,3,-1,1)
		nevals=collect(linspace(-1,1,65))
		basis1 = full(getBasis(nevals,version1))
		y=f3(nevals)
		coef1 = basis1 \ y
		base1=full(getBasis(points,version1))
		estimation1=base1 * coef1
		knots=zeros(13)
		if n%2==0
			#we want to use linspace with a pair number of points to not have any zero.
			non0knots=collect(linspace(-1,1,13-n-1))
		else
			non0knots=collect(linspace(-1,1,13-n))
		end
		for i in 1:length(non0knots)
			knots[i]=non0knots[i]
		end
		knots=sort(knots)
		#And same as before, we evaluate the coefficient at the evaluation nodes and compute an estimation for 500 points
		params = SplineParams(knots,0,3)
		basis2=CompEcon.evalbase(params,nevals)[1]
		coef2 = basis2 \ y
		base2=CompEcon.evalbase(params,points)[1]
		estimation2=base2 * coef2
		#and plotting!
		figure("Question 5")
		subplot(311)
		plot(points,truth,color="blue")
		title("true values")
		subplot(312)
		plot(points, estimation1, color="red")
		plot(points, estimation2, color="green")
		title("Estimated values: red no or 1 zero knots, green $n zero knots")
		subplot(313)
		plot(points, truth .- estimation1, color="red")
		plot(points, truth .- estimation2, color="green")
		title("Deviation of the estimations")
	end

###################################END OF QUESTION 5###########################


		# function to run all questions
	function runall()
		println("running all questions of HW-funcapprox:")
		println("Question 1:")
		fapprox=q1(15)
		x=linspace(-3,3,50)
		y_approx=[]
		y_exact=zeros(50)
		for i in 1:50
			#our interpolator returns an array
			y_approx=push!(y_approx,fapprox(x[i]))
			y_exact[i]=f(x[i])
		end
		deviation=y_exact .- y_approx
		figure("Question 1")
		plot(x,y_approx,color="blue")
		plot(x,y_exact, color="red")
		plot(x,deviation, color="green")
		title("The true function (red), the approximated one (blue) and the difference (green)")
		xlabel("x")
		ylabel("y")
##############################QUESTION 2########################################
		println("Question 2")
		fapprox=q2(15)
		x=linspace(-3,3,50)
		y_approx=[]
		y_exact=zeros(50)
		for i in 1:50
			#our interpolator returns an array
			y_approx=push!(y_approx,fapprox(x[i]))
			y_exact[i]=f(x[i])
		end
		deviation=y_exact .- y_approx
		figure("Question 2")
		plot(x,y_approx,color="blue")
		plot(x,y_exact, color="red")
		plot(x,deviation, color="green")
		title("The true function (red), the approximated one (blue) and the difference (green)")
		xlabel("x")
		ylabel("y")
#############################QUESTION 3########################################
		q3()
#############################QUESTION 4a#######################################
		q4a()
#############################QUESTION 4b#######################################
		q4b()
#############################QUESTION 5########################################
		q5()
	end
end
