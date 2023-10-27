using System;

// Recreation of the UnityEngine.Vector3, so it can be used in other thread
public struct Vector4d
{

	public double x, y, z, w;

	public static readonly Vector4d zero = new Vector4d(0, 0, 0, 0);
	public static readonly Vector4d one = new Vector4d(1, 1, 1, 1);

	public Vector4d(double x, double y)
	{
		this.x = x;
		this.y = y;
		this.z = 0;
		this.w = 1;
	}
	public Vector4d(double x, double y, double z, double w)
	{
		this.x = x;
		this.y = y;
		this.z = z;
		this.w = w;
	}

	public double this[int i]
	{
		get
		{
			double res = 0;
			if (i == 0)
			{
				res = this.x;
			}
			else if (i == 1)
			{
				res = this.y;
			}
			else if (i == 2)
			{
				res = this.z;
			}
			else if(i == 3)
			{
				res = this.w;
			}
			return res;
		}
		set
		{
			if (i == 0)
			{
				this.x = (double)value;
			}
			else if (i == 1)
			{
				this.y = (double)value;
			}
			else if (i == 2)
			{
				this.z = (double)value;
			}
			else if (i == 3)
			{
				this.w = (double)value;
			}
		}
	}

	public double magnitude
	{
		get
		{
			return (double)Math.Sqrt(x * x + y * y + z * z + w * w);
		}
	}

	public void Normalize()
	{
		double magnitude = this.magnitude;
		x /= magnitude;
		y /= magnitude;
		z /= magnitude;
		w /= magnitude;
	}

	public static Vector4d Normalize(Vector4d a)
	{
		double magnitude = a.magnitude;
		return new Vector4d(a.x / magnitude, a.y / magnitude, a.z / magnitude, a.w / magnitude);
	}

	public override bool Equals(object other)
	{
		if (!(other is Vector4d))
		{
			return false;
		}
		Vector4d v = (Vector4d)other;
		return x == v.x && y == v.y && z == v.z && w == v.w;
	}

	public override string ToString()
	{
		return string.Format("[Vector4d]" + x + "," + y + "," + z + "," + w);
	}

	public override int GetHashCode()
	{
		return x.GetHashCode() ^ y.GetHashCode() << 2;
	}

	public double DistanceSquare(Vector4d v)
	{
		return Vector4d.DistanceSquare(this, v);
	}
	public static double DistanceSquare(Vector4d a, Vector4d b)
	{
		double cx = b.x - a.x;
		double cy = b.y - a.y;
		double cz = b.z - a.z;
		double cw = b.w - a.w;
		return cx * cx + cy * cy + cz * cz + cw * cw;
	}

	public static bool operator ==(Vector4d a, Vector4d b)
	{
		return a.x == b.x && a.y == b.y && a.z == b.z  && a.w == b.w;
	}

	public static bool operator !=(Vector4d a, Vector4d b)
	{
		return a.x == b.x || a.y == b.y || a.z == b.z || a.w == b.w;
	}

	public static Vector4d operator -(Vector4d a, Vector4d b)
	{
		return new Vector4d(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
	}

	public static Vector4d operator +(Vector4d a, Vector4d b)
	{
		return new Vector4d(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
	}

	public static Vector4d operator *(Vector4d a, int i)
	{
		return new Vector4d(a.x * i, a.y * i, a.z * i, a.w * i);
	}

	public static Vector4d operator *(Vector4d a, double i)
	{
		return new Vector4d(a.x * i, a.y * i, a.z * i, a.w * i);
	}

	//public static Vector4d Min(Vector4d a, Vector4d b)
	//{
	//	return new Vector4d(Math.Min(a.x, b.x), Math.Min(a.y, b.y));
	//}
	//public static Vector2f Max(Vector2f a, Vector2f b)
	//{
	//	return new Vector2f(Math.Max(a.x, b.x), Math.Max(a.y, b.y));
	//}
}
