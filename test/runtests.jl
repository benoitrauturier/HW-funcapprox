module Test
  using funcapp
  using FactCheck
  using PyPlot
##############################TEST Q1 #########################################
  facts("testing deviation Q1") do
    fapprox=funcapp.q1(15)
    x=linspace(-3.0,3.0,500)
    maxim=0.0
    for i in 1:500
      val=abs(fapprox(x[i])-funcapp.f(x[i]))
      if val[1] > maxim
        maxim = val[1]
      end
    end
    println(maxim)
    @fact float(maxim)<1e-8 --> true
  end
################################TEST Q2#########################################
  facts("testing deviation Q2") do
    fapprox=funcapp.q2(15)
    x=linspace(-3,3,500)
    maxim=0.0
    for i in 1:500
      val=abs(fapprox(x[i])-funcapp.f(x[i]))
      if val[1] > maxim
        maxim = val[1]
      end
    end
    println(maxim)
    @fact float(maxim)<1e-9 --> true
  end
end
