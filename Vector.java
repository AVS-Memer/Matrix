import java.util.Arrays;
import java.util.function.*;
public class Vector {
  public static final Vector EMPTY = new Vector();
  double[] values; int size;
  public Vector() {
    values = new double[0];
    size = 0;
  }
  public Vector(int i) {
    if (i < 0) {
      throw new IllegalArgumentException("Invalid Vector Dimensions");
    }
    size = i;
    values = new double[i];
    Arrays.fill(values, 0);
  }
  public Vector(double[] v) {
    size = v.length;
    values = Arrays.copyOf(v, size);
  }
  public static Vector fill(int a, double x) {
    return new Vector(a).map(m -> x);
  }
  public static Vector zeros(int a) {
    return new Vector(a);
  }
  public static Vector ones(int a) {
    return new Vector(a).map(m -> 1.0);
  }
  public static Vector random(int a) {
    return new Vector(a).map(m -> Math.random());
  }
  public static Vector fromMatrix(Matrix n) {
    if (n.isEmpty()) return new Vector();
    else if (n.size[0] == 1) {
      return new Vector(n.values[0]);
    } else if (n.size[1] == 1) {
      Vector v = new Vector(n.size[0]);
      for (int i = 0; i < n.size[0]; i++) v.values[i] = n.values[i][0];
      return v;
    } else {
      throw new IllegalArgumentException("Invalid Matrix Dimensions for Vector Creation");
    }
  }
  public int size() {
    return size;
  }
  public Vector augment(double m2) {
    double[] r = new double[size+1];
    System.arraycopy(values, 0, r, 0, size);
    r[size] = m2;
    return new Vector(r);
  }
  public Vector augment(Vector m2) {
    double[] r = new double[size+m2.size];
    System.arraycopy(values, 0, r, 0, size);
    System.arraycopy(m2.values, 0, r, size, m2.size);
    return new Vector(r);
  }
  public double sum() {
    double sum = 0;
    for (int i = 0; i < size; i++) sum += values[i];
    return sum;
  }
  public Vector submatrix(int start, int end) {
    if (start < 0 || end > size || start > end)
      throw new IllegalArgumentException("Invalid subvector bounds");
    double[] sub = new double[end-start];
    System.arraycopy(values, start, sub, 0, end-start);
    return new Vector(sub);
  }
  public Vector submatrix(int end) {
    return submatrix(0, end);
  }
  public double get(int r) {
    if (r >= size) throw new IllegalArgumentException("Invalid index");
    return values[r];
  }
  public void set(int r, double x) {
    if (r >= size) throw new IllegalArgumentException("Invalid row index");
    values[r] = x;
  }
  public void swapElements(int r1, int r2) {
    if (r1 >= size || r2 >= size) throw new IllegalArgumentException("Invalid indexes for swapping");
    if (r1 != r2) {
      double temp = values[r1];
      values[r1] = values[r2];
      values[r2] = temp;
    }
  }
  public double mean() {
    return this.sum()/(size);
  }
  public double median() {
    double[] allValues = Arrays.copyOf(values,size);
    Arrays.sort(allValues);
    return (allValues.length % 2 == 1) ? allValues[allValues.length / 2] : (allValues[allValues.length / 2 - 1] + allValues[allValues.length / 2]) / 2.0;
  }
  public Vector center() {
    return this.subtract(this.mean());
  }
  public double std() {
    return Math.sqrt(this.center().hadamard(this.center()).mean());
  }
  public boolean isEqual(Vector m2) {
    if (size != m2.size) return false;
    for (int i = 0; i < size; i++) {
      if (Math.abs(values[i] - m2.values[i]) > 1e-9) return false;
    }
    return true;
  }
  public boolean isEmpty() {
    return this.isEqual(Vector.EMPTY);
  }
  public Vector add(double m2) {
    return this.map(a -> a+m2);
  }
  public Vector add(Vector m2) {
    if (size != m2.size) {
      throw new IllegalArgumentException("Vectors must be same size for addition.");
    }
    return this.map(m2,(a,b) -> a+b);
  }
  public Vector subtract(double m2) {
    return this.map(a -> a-m2);
  }
  public Vector subtract(Vector m2) {
    if (size != m2.size) {
      throw new IllegalArgumentException("Vectors must be same size for subtraction.");
    }
    return this.map(m2,(a,b) -> a-b);
  }
  public Vector multiply(double number) {
    return this.map(a -> a*number);
  }
  public Vector negate() {
    return this.map(a -> -a);
  }
  public Vector hadamard(Vector m2) {
    return this.map(m2, (a,b) -> a*b);
  }
  public Vector multiply(Matrix m2) {
    if (m2.size[0] != size) {
      throw new IllegalArgumentException("Matrix columns must match vector length.");
    }
    double[] v = new double[m2.size[1]];
    for (int i = 0; i < m2.size[1]; i++) {
      v[i] = 0;
      for (int j = 0; j < m2.size[0]; j++) {
        v[i] += m2.values[j][i] * values[j];
      }
    }
    return new Vector(v);
  }
  public Vector divide(double number) {
    return this.map(a -> a/number);
  }
  public double norm() {
    return Math.sqrt(this.dotProduct(this));
  }
  public Vector normalize() {
    return this.divide(this.norm());
  }
  public double distance(Vector m2) {
    return this.subtract(m2).norm();
  }
  public double dotProduct(Vector m2) {
    if (size != m2.size) {
      throw new IllegalArgumentException("Vectors must be same size for dot product.");
    }
    double sum = 0;
    for (int i = 0; i < size; i++) sum += values[i]*m2.values[i];
    return sum;
  }
  public double angleBetweenVectors(Vector m2) {
    return Math.acos(this.normalize().dotProduct(m2.normalize()));
  }
  public Vector crossProduct(Vector m2) {
    if (size != 3 || m2.size != 3) {
      throw new IllegalArgumentException("Cross product is only defined for 3D vectors.");
    }
    return new Vector(new double[]{this.values[1]*m2.values[2]-this.values[2]*m2.values[1], this.values[2]*m2.values[0]-this.values[0]*m2.values[2], this.values[0]*m2.values[1]-this.values[1]*m2.values[0]});
  }
  public double tripleScalarProduct(Vector m2, Vector m3) {
    return this.dotProduct(m2.crossProduct(m3));
  }
  public Vector tripleVectorProduct(Vector m2, Vector m3) {
    return this.crossProduct(m2.crossProduct(m3));
  }
  public double linearCombination(double[] m2) {
    return this.hadamard(new Vector(m2)).sum();
  }
  public double linearCombination(Vector m2) {
    return this.hadamard(m2).sum();
  }
  public Vector projection(Vector m2) {
    return m2.multiply(this.dotProduct(m2)/m2.dotProduct(m2));
  }
  public Vector reflection(Vector m2) {
    return this.projection(m2).multiply(2).subtract(this);
  }
  public Vector parallelComponent(Vector m2) {
    return this.projection(m2);
  }
  public Vector perpendicularComponent(Vector m2) {
    return this.subtract(this.projection(m2));
  }
  public Matrix outerProduct(Vector m2) {
    double[][] v = new double[size][m2.size];
    for (int i = 0; i < size; i++) {
      System.arraycopy(m2.multiply(values[i]).values,0,v[i],0,m2.size);
    }
    return new Matrix(v);
  }
  public Vector kroneckerProduct(Vector m2) {
    double[] v = new double[size*m2.size];
    for (int i = 0; i < size; i++) {
      System.arraycopy(m2.multiply(values[i]).values,0,v,i*m2.size,m2.size);
    }
    return new Vector(v);
  }
  public Vector midpoint(Vector m2) {
    return this.map(m2,(a,b) -> (a+b)/2);
  }
  public double max() {
    if (size == 0) throw new IllegalArgumentException("Empty vector does not contain maximum value");
    double m = values[0];
    for (double v : values) m = Math.max(v,m);
    return m;
  }
  public double min() {
    if (size == 0) throw new IllegalArgumentException("Empty vector does not contain minimum value");
    double m = values[0];
    for (double v : values) m = Math.min(v,m);
    return m;
  }
  public Vector softmax() {
    double[] expVals = new double[size];
    double max = this.max();
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
      expVals[i] = Math.exp(values[i] - max);
      sum += Math.exp(values[i] - max);
    }
    return new Vector(expVals).divide(sum);
  }
  public Vector conditionalMap(Predicate<Double> condition, Function<Double, Double> ifTrue, Function<Double, Double> ifFalse) {
    return this.map(a -> condition.test(a) ? ifTrue.apply(a) : ifFalse.apply(a));
  }
  public Vector map(Function<Double, Double> f) {
    double[] r = new double[size];
    for (int i = 0; i < size; i++) r[i] = f.apply(values[i]);
    return new Vector(r);
  }
  public Vector map(Vector m2, BiFunction<Double, Double, Double> func) {
    if (size != m2.size) {
      throw new IllegalArgumentException("Vectors must be the same size for elementwise map.");
    }
    double[] r = new double[size];
    for (int i = 0; i < size; i++) r[i] = func.apply(values[i], m2.values[i]);
    return new Vector(r);
  }
  public Vector applyActivation(ActivationFunction f) {
    return this.map(f.get());
  }
  @Override
  public Vector clone() {
    return new Vector(Arrays.copyOf(values, size));
  }
  @Override
  public String toString() {
    return this.toString(0);
  }
  public String toString(int decimals) {
    int maxLength = 0;
    for (int i = 0; i < size; i++) maxLength = Math.max(maxLength, String.format("%."+decimals+"f", values[i]).length());
    StringBuilder s = new StringBuilder();
    s.append("< ");
    for (int i = 0; i < size; i++) {
      s.append(String.format("%"+maxLength+"."+decimals+"f", values[i]));
      if (i < size - 1) s.append("  ");
    }
    s.append(" >");
    return s.toString();
  }
}
