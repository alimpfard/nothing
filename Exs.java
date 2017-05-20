package net.anothertest;
import java.lang.reflect.*;
import java.util.*;
import java.util.stream.*;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.image.BufferedImage;


public class Exs{
	interface Q{
		public void run(String[] args);
		}
	static class Q1 implements Q{
		class Rational{
			protected int numerator, denumerator;
			public Rational(int num, int den){numerator=num; denumerator=den;}
			private Rational(){}
			private Rational negate(){numerator*=-numerator; return this;}
			private Rational invert(){numerator^=denumerator; denumerator^=numerator; numerator^=denumerator; return this;}
			private Rational copy(){
				Rational r=new Rational();
				r.numerator=this.numerator;
				r.denumerator=this.denumerator;
				return r;
			}
			public Rational plus(Rational b){
				Rational s = b.copy();
				s.denumerator *= denumerator;
				s.numerator = s.numerator*denumerator + numerator*s.denumerator;
				return s;
			}
			public Rational minus(Rational b){
				Rational s = b.copy().negate();
				s.denumerator *= denumerator;
				s.numerator = s.numerator*denumerator + numerator*s.denumerator;
				return s;
			}
			public Rational times(Rational b){
				Rational r = b.copy();
				r.denumerator*=denumerator; r.numerator*=denumerator;
				return r;
			}
			public Rational divides(Rational b){
				Rational r = b.copy().invert();
				r.denumerator*=denumerator; r.numerator*=denumerator;
				return r;
			}
			public String toString(){
				return String.format("(%s/%s)",numerator,denumerator);
			}
		}
		@Override
		public void run(String[] args){
			Rational a = new Rational(2,3);
			Rational b = new Rational(3,4);
			System.out.println(a.times(b).plus(a).divides(a).minus(b));
		}
	}
	static class Q2 implements Q{
		strictfp class Interval{
			protected double left, right;
			protected Boolean empty = false;
			public Interval(double l, double r){left=l; right=r; empty=l>r;}
			private Interval(){}
			public Boolean contains(double x){return !empty&&x<=right&&x>=left;}
			public Boolean intersects(Interval b){
				Boolean intersect=false;
				intersect = contains(b.left)||contains(b.right);
				return intersect&&!empty;
			}
			public String toString(){
				return (!empty?String.format("[%s,%s]",left,right):"Empty Interval");
			}
		}
		@Override
		public void run(String[] args){
			Scanner s = new Scanner(System.in);
			int mkN = 0;
			double filter = Double.parseDouble(args[0]);
			List<Interval> intervals = new LinkedList<Interval>();
			double x[] = new double[2];
			while (s.hasNextDouble()){
				x[mkN++]=s.nextDouble();
				if(mkN==1){
					intervals.add(new Interval(x[0], x[1]));
					mkN=0;
				}
			}
			intervals.stream().filter(i->i.contains(filter)).forEach(i->System.out.println(i));
		}
	}
	static class Q3 implements Q{
		class Point{
			protected final double x, y;
			public Point(double _x, double _y){
				x = _x; y = _y;
			}
			public double distanceTo(Point p){
				return Math.sqrt((p.x-x)*(p.x-x)  + (p.y-y)*(p.y-y));
			}
			public String toString(){
				return String.format("(%s,%s)",x,y);
			}
		}
		@Override
		public void run(String[] args){
			Point a = new Point(0,1), b = new Point(2,4.32d);
			System.out.println(String.format("The distance between points %s and %s is %s",a,b,a.distanceTo(b)));
		}
	}
	static class Q4 implements Q{
		static class Complex {
    protected final double re;
    protected final double im;
    public Complex(double real, double imag) { re = real; im = imag; }
		public Complex copy(){return new Complex(re,im);}
    public String toString() {
        if (im == 0) return re + "";
        if (re == 0) return im + "i";
        if (im <  0) return re + " - " + (-im) + "i";
        return re + " + " + im + "i";
    }
    public double abs() {
        return Math.hypot(re, im);
    }
    public Complex plus(Complex b) {
        Complex a = this;
        double real = a.re + b.re;
        double imag = a.im + b.im;
        return new Complex(real, imag);
    }
		public Complex minus(Complex b) {
        Complex a = this;
        double real = a.re - b.re;
        double imag = a.im - b.im;
        return new Complex(real, imag);
    }
    public Complex times(Complex b) {
        Complex a = this;
        double real = a.re * b.re - a.im * b.im;
        double imag = a.re * b.im + a.im * b.re;
        return new Complex(real, imag);
    }
    public Complex conjugate() {
        return new Complex(re, -im);
    }
		public Complex reciprocal() {
		double scale = re*re + im*im;
		return new Complex(re / scale, -im / scale);
		}
    public double re() { return re; }
    public double im() { return im; }
    public Complex divides(Complex b) {
        Complex a = this;
        return a.times(b.reciprocal());
    }
    public Complex exp() {
        return new Complex(Math.exp(re) * Math.cos(im), Math.exp(re) * Math.sin(im));
    }
		public double phase() {
			 return Math.atan2(im, re);
		 }
		 public Collection<Complex> nthRoot(int n) throws IllegalArgumentException {
	 if (n <= 0) {
		 throw new IllegalArgumentException("The value for the nth root has to be positive!");
	 }
		 Collection<Complex> result = new ArrayList<Complex>();
	 // nth root of abs
	 double nthRootOfAbs = Math.pow( abs() , 1.0/n );
	 // Compute nth roots of complex number with k=0,1,...n-1
		 for (int k=0; k<n;k++) {
			 // inner part
			 double innerPart = (phase() + k * 2 * Math.PI) / n;
			 double realPart = nthRootOfAbs *  Math.cos ( innerPart );
			 double imaginaryPart = nthRootOfAbs *  Math.sin ( innerPart );
			 result.add(new Complex(realPart, imaginaryPart));
		 }
		 return result;
 }
		//Who gives a rat's ass?
		public Complex pow(int b){
			Complex c = this.copy();
			final Complex c0 = c.copy();
			for(int i=0; i<b-1; i++)
			c = c.times(c0);
			return c;
		}

	}
	@Override
	public void run(String[] args){
		Complex a = new Complex(0,10), b = new Complex(10,3.5d);
		System.out.println(String.format("%s, %s, to the third power: %s, %s",a,b,a.pow(3), b.pow(3)));
	}
}
	static class Q5 implements Q{
		class RootOfUnity extends Q4.Complex{
			public RootOfUnity(double x, double y){super(x,y);}
			public Collection<Q4.Complex> solve(int n){
				return nthRoot(n);
			}
		}
		public void run(String[] args){
			RootOfUnity rou = new RootOfUnity(Double.parseDouble(args[0]), Double.parseDouble(args[1]));
			String res = rou.solve(Integer.parseInt(args[2])).stream().map(i->i+"").collect(Collectors.joining(",\n"));
			System.out.println(String.format("roots for %s:\n%s",rou,res));
		}
	}
	static class Q6 implements Q {
		class Complex extends Q4.Complex{
			public Complex(double x, double y){super(x,y);}
		}
		@Override
		public void run(String[] args){
			Complex x = new Complex(0,3), y = new Complex(3,4);
			System.out.println(x.pow(2).divides(y).minus(y).conjugate().phase());
		}
	}
	static class Q7 implements Q{
		@Override
		public void run(String[] args){
			System.out.println("You have defined a void method, named Complex--which won't compile anyway--,\n if you wanted a %tor, you should've ommited the `void` type.");
		}
	}
	static class Q8 implements Q{
		class Charge{
			protected double rx, ry, q;
			public Charge(double rx, double ry, double q){this.rx=rx; this.ry=ry; this.q=q;}
			public double potentialAt(double x, double y){
				double k = 8.99e09;
				double dx = x - rx, dy = y - ry;
				return k*q/Math.sqrt(dx*dx+dy*dy);
			}
			public String toString(){
				return String.format("(charge%s at %s,%s)",q,rx,ry);
			}
		}
		public Color getColor(double v){
			v = 128 + v / 2.0e10;
			int t = 0;
			if      (v <   0) t = 0;
			else if (v > 255) t = 255;
			else              t = (int) v;
			int gray = (t * 37) % 255;
			return new Color(gray,gray,gray);
		}
		@Override
		public void run(String[] args){
			ArrayList<Charge> a = new ArrayList<Charge>();
			Scanner s = new Scanner(System.in);
			while (s.hasNextDouble()){
				double x=s.nextDouble(), y=s.nextDouble(), z=s.nextDouble();
				a.add(new Charge(x,y,z));
			}
			int SIZE = 512;
			Picture pic = new Picture(SIZE, SIZE);
			for (int i = 0; i < SIZE; i++) {
			for (int j = 0; j < SIZE; j++) {
			double V = 0.0;
			for (int k = 0; k < a.size(); k++) {
			double x = 1.0 * i / SIZE;
			double y = 1.0 * j / SIZE;
			V += a.get(k).potentialAt(x, y);
			}
			Color color = getColor(V);
			pic.set(i, SIZE-1-j, color);
			}
			}
			pic.show();
		}
	}
	static class Q9 implements Q{
		public class Quaternion { //SHAMELESSLY STOLEN FROM UNITY
    protected double x0, x1, x2, x3;
    public Quaternion(double _x0, double _x1, double _x2, double _x3) {
        x0 = _x0;
        x1 = _x1;
        x2 = _x2;
        x3 = _x3;
    }
    public String toString() {
        return String.format("(Quat%s + %si + %sj + %sk)",x0,x1,x2,x3);
    }
    public double norm() {
        return Math.sqrt(x0*x0 + x1*x1 +x2*x2 + x3*x3);
    }

    public Quaternion conjugate() {
        return new Quaternion(x0, -x1, -x2, -x3);
    }
    public Quaternion plus(Quaternion b) {
        return new Quaternion(x0+b.x0, x1+b.x1, x2+b.x2, x3+b.x3);
    }
    public Quaternion times(Quaternion b) {
        double y0 = x0*b.x0 - x1*b.x1 - x2*b.x2 - x3*b.x3,
        			 y1 = x0*b.x1 + x1*b.x0 + x2*b.x3 - x3*b.x2,
        			 y2 = x0*b.x2 - x1*b.x3 + x2*b.x0 + x3*b.x1,
        			 y3 = x0*b.x3 + x1*b.x2 - x2*b.x1 + x3*b.x0;
        return new Quaternion(y0, y1, y2, y3);
    }
    public Quaternion inverse() {
        double d = x0*x0 + x1*x1 + x2*x2 + x3*x3;
        return new Quaternion(x0/d, -x1/d, -x2/d, -x3/d);
    }
    public Quaternion divides(Quaternion b) {
        return this.times(b.inverse());
    }
	}
	@Override
	public void run(String[] args){
		Quaternion q = new Quaternion(0,1,2,3);
		System.out.println(q.times(q).plus(q).inverse().divides(q).conjugate().times(q));
	}
}
	static class Q10 implements Q{
		static class Turtle{
			protected double x, y; // turtle is at (x, y)
			protected double angle; // facing this direction
			Random r = new Random(); //with this random gen -_-
			public Turtle(double x0, double y0, double a0) { x = x0; y = y0; angle = a0; }
			public void left(double delta) {angle += delta;}
			public void right(double d){angle-=d;}
			public void forward(double d) {
			double oldx = x;
			double oldy = y;
			x += d * Math.cos(Math.toRadians(angle));
			y += d * Math.sin(Math.toRadians(angle));
			StdDraw.line(oldx, oldy, x, y);
			}
			public int choice(int[] from){
				return from[r.nextInt(from.length)];
			}
		}
		// dragon curve of order n
		public void dragon(int n, Turtle turtle) {
				if (n == 0) {
						turtle.forward(1.0);
				}

				else {
						dragon(n-1,turtle);
						turtle.left(90);
						revdrag(n-1,turtle);
				}
		}
		// reverse dragon curve of order n
		public void revdrag(int n, Turtle turtle) {
				if (n == 0) {
						turtle.forward(1.0);
				}
				else {
						dragon(n-1,turtle);
						turtle.left(-90);
						revdrag(n-1,turtle);
				}
			}
		@Override
		public void run(String[] args){ //Some Magic Numbers here... They say these are here to figure out the size...i ain't no Math-Geek though.
			int[] left  = { 0, 0, 0, 2, 4, 5, 5,  5,  5,  5, 10, 42, 74, 81,  85,  85 };
      int[] right = { 1, 1, 1, 1, 1, 1, 2, 10, 18, 21, 21, 21, 21, 21,  57, 170 };
      int[] up    = { 0, 1, 2, 2, 2, 2, 2,  2,  5, 21, 37, 42, 42, 42,  42,  42 };
      int[] down  = { 0, 0, 0, 0, 1, 5, 9, 10, 10, 10, 10, 10, 23, 85, 149, 165 };
			Turtle t = new Turtle(0,0,0);
			int order = Integer.parseInt(args[0]);
			order = order>15?15:order;
			double size = Math.max(left[order]+right[order],up[order]+down[order]),
			       x = (right[order]-left[order])/2,
			       y = (up[order]-down[order])/2;
		 StdDraw.setXscale(x - size/2, x + size/2);
     StdDraw.setYscale(y - size/2, y + size/2);
		 dragon(order,t);
		}
	}
	static class Q11 implements Q{
		class Turtle extends Q10.Turtle{public Turtle(int x, int y, int z){super(x,y,z);}}
			// Hilbert curve
	    private void hilbert(int n, Turtle turtle) {
	        if (n == 0) return;
	        turtle.left(90);
	        treblih(n-1,turtle);
	        turtle.forward(1.0);
	        turtle.left(-90);
	        hilbert(n-1,turtle);
	        turtle.forward(1.0);
	        hilbert(n-1,turtle);
	        turtle.left(-90);
	        turtle.forward(1.0);
	        treblih(n-1,turtle);
	        turtle.left(90);
	    }

	    // reverse Hilbert
	    public void treblih(int n, Turtle turtle) {
	        if (n == 0) return;
	        turtle.left(-90);
	        hilbert(n-1,turtle);
	        turtle.forward(1.0);
	        turtle.left(90);
	        treblih(n-1,turtle);
	        turtle.forward(1.0);
	        treblih(n-1,turtle);
	        turtle.left(90);
	        turtle.forward(1.0);
	        hilbert(n-1,turtle);
	        turtle.left(-90);
	    }
		@Override
		public void run(String[] args){
			Turtle t = new Turtle(0,0,0);
			int n = Integer.parseInt(args[0]);
			double max = Math.pow(2, n);
      StdDraw.setXscale(0, max);
      StdDraw.setYscale(0, max);
			hilbert(n,t);
		}
	}
	static class Q12 implements Q{
		private final double Angle = 60;
		class Turtle extends Q10.Turtle{
			protected double size;
			public Turtle(double x, double y, double z){super(x,y,z);}
			public void forward(){forward(size);}
			public void setSize(double s){size=s;}
		}
		private void gosper(int n, Turtle turtle) {
			if (n == 0) turtle.forward();
			else {
					turtle.right(19.106605350869094394);
					gosper(n-1,turtle);
					turtle.left(Angle);
					gosper(n-1,turtle);
					turtle.right(Angle);
					gosper(n-1,turtle);
				  turtle.left(19.106605350869094394);
			}
		}
		private void gisland(int n, Turtle turtle){
			turtle.setSize((7.0 / 16.0) / Math.pow(Math.sqrt(7.0), n));
			gosper(n,turtle);
        turtle.left(60);
        gosper(n,turtle);
        turtle.left(60);
        gosper(n,turtle);
        turtle.left(60);
        gosper(n,turtle);
        turtle.left(60);
        gosper(n,turtle);
        turtle.left(60);
        gosper(n,turtle);
        turtle.left(60);
		}
	@Override
	public void run(String[] args){
		Turtle t = new Turtle(0,0, 0.0);
		int n = Integer.parseInt(args[0]);
		StdDraw.setScale(-10/n,10/n);
		gisland(n,t);
	}
	}
	static class Q13 implements Q{
		class Complex extends Q4.Complex{
			public Complex(double x, double y){super(x,y);}
		}
		public  Color newton(Q4.Complex z) {
				double EPSILON = 0.000001;
				Complex four = new Complex(4, 0);
				Complex one  = new Complex(1, 0);
				Complex root1 = new Complex(1, 0);
				Complex root2 = new Complex(-1, 0);
				Complex root3 = new Complex(0, 1);
				Complex root4 = new Complex(0, -1);
				for (int i = 0; i < 100; i++) {
						Q4.Complex f  = z.pow(4).minus(one);
						Q4.Complex fp = four.times(z.pow(3));
						z = z.minus(f.divides(fp));
						if (z.minus(root1).abs() <= EPSILON) return Color.WHITE;
						if (z.minus(root2).abs() <= EPSILON) return Color.RED;
						if (z.minus(root3).abs() <= EPSILON) return Color.GREEN;
						if (z.minus(root4).abs() <= EPSILON) return Color.BLUE;
				}
				return Color.BLACK;
			}
			@Override
			public void run(String[] args){
				int n = Integer.parseInt(args[0]);
				double xmin   = -1.0,
				 			 ymin   = -1.0,
				 			 width  =  2.0,
				 			 height =  2.0;
				Picture picture = new Picture(n, n);
				for (int col = 0; col < n; col++) {
						for (int row = 0; row < n; row++) {
								double x = xmin + col * width  / n,
								 			 y = ymin + row * height / n;
								Complex z = new Complex(x, y);
								Color color = newton(z);
								picture.set(col, row, color);
						}
				}
				picture.show();
				}
			}
	static class Q14 implements Q{
		public int mandel(Q4.Complex z0, int d) {
        	Q4.Complex z = z0;
        	for (int t = 0; t < d; t++) {
        	    if (z.abs() > 2.0) return t;
        	    z = z.pow(2).plus(z0);
        	}
        	return d;
        	}
        @Override
        public void run(String[] args){
        	Scanner s = new Scanner(System.in);
        	double xc = Double.parseDouble(args[0]),
        	       yc = Double.parseDouble(args[1]),
        	       size = Double.parseDouble(args[2]);
        	int n = 512,
        	ITERS = 256;
        	// read color map
        	Color[] colors = new Color[ITERS];
        	for (int t = 0; t < ITERS; t++) {
        		int r = s.nextInt(),
        		    g = s.nextInt(),
        		    b = s.nextInt();
        		colors[t] = new Color(r, g, b); }
        	// Mandelbrot set
        	Picture picture = new Picture(n, n);
        	for (int col = 0; col < n; col++) {
        		for (int row = 0; row < n; row++) {
        			double x = xc - size/2 + size*col/n,
        			       y = yc - size/2 + size*row/n;
        			Q4.Complex z0 = new Q4.Complex(x, y);
        			int t = mandel(z0, ITERS - 1);
        			picture.set(col, n-1-row, colors[t]);
        			}
        		}
        	picture.show();
        	}
        }
	static class Q15 implements Q{
		public int julia(Q4.Complex c, Q4.Complex z, int maximumIterations) {
        	for (int t = 0; t < maximumIterations; t++) {
            	if (z.abs() > 2.0) return t;
            	z = z.times(z).plus(c);
        	}
        	return maximumIterations - 1;
    }
    @Override
    public void run(String[] args){
        Scanner s = new Scanner(System.in);
        double real = Double.parseDouble(args[0]);      // such Math much Wow
        double imag = Double.parseDouble(args[1]);
        Q4.Complex c = new Q4.Complex(real, imag);
        double xmin   = -2.0,
               ymin   = -2.0,
               width  =  4.0,
               height =  4.0;

        int n = 512;
        int ITERS  = 256;

        Color[] colors = new Color[ITERS];
        for (int t = 0; t < ITERS; t++) {
            int r = s.nextInt();
            int g = s.nextInt();
            int b = s.nextInt();
            colors[t] = new Color(r, g, b);
        }
        Picture picture = new Picture(n, n);

        for (int col = 0; col < n; col++) {
            for (int row = 0; row < n; row++) {
                double x = xmin + col * width / n,
                       y = ymin + row * height / n;
                Q4.Complex z = new Q4.Complex(x, y);
                int t = julia(c, z, ITERS);
                picture.set(col, row, colors[t]);
            }
        }
        picture.show();
        }

	} 
	

	public Exs(){}
	public static void main(String args[]){
	Q question=null;
	try{Class c = Class.forName("net.anothertest.Exs$Q"+args[0]);
	question = (Q)(c.newInstance());}
	catch(Exception e){System.out.println(e);}
	question.run(Arrays.asList(args).subList(1,args.length).toArray(new String[0]));
	}

}
