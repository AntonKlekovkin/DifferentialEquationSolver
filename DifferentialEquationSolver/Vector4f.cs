using System;

// Recreation of the UnityEngine.Vector3, so it can be used in other thread
public struct Vector4f
{

	public float x, y, z, w;

	public static readonly Vector4f zero = new Vector4f(0, 0, 0, 0);
	public static readonly Vector4f one = new Vector4f(1, 1, 1, 1);

	public Vector4f(float x, float y)
	{
		this.x = x;
		this.y = y;
		this.z = 0;
		this.w = 1;
	}
	public Vector4f(float x, float y, float z, float w)
	{
		this.x = x;
		this.y = y;
		this.z = z;
		this.w = w;
	}
	public Vector4f(double x, double y, double z, double w)
	{
		this.x = (float)x;
		this.y = (float)y;
		this.z = (float)z;
		this.w = (float)w;
	}

	public double this[int i]
	{
		get
		{
			float res = 0;
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
				this.x = (float)value;
			}
			else if (i == 1)
			{
				this.y = (float)value;
			}
			else if (i == 2)
			{
				this.z = (float)value;
			}
			else if (i == 3)
			{
				this.w = (float)value;
			}
		}
	}

	public float magnitude
	{
		get
		{
			return (float)Math.Sqrt(x * x + y * y + z * z + w * w);
		}
	}

	public void Normalize()
	{
		float magnitude = this.magnitude;
		x /= magnitude;
		y /= magnitude;
		z /= magnitude;
		w /= magnitude;
	}

	public static Vector4f Normalize(Vector4f a)
	{
		float magnitude = a.magnitude;
		return new Vector4f(a.x / magnitude, a.y / magnitude, a.z / magnitude, a.w / magnitude);
	}

	public override bool Equals(object other)
	{
		if (!(other is Vector4f))
		{
			return false;
		}
		Vector4f v = (Vector4f)other;
		return x == v.x && y == v.y && z == v.z && w == v.w;
	}

	public override string ToString()
	{
		return string.Format("[Vector4f]" + x + "," + y + "," + z + "," + w);
	}

	public override int GetHashCode()
	{
		return x.GetHashCode() ^ y.GetHashCode() << 2;
	}

	public float DistanceSquare(Vector4f v)
	{
		return Vector4f.DistanceSquare(this, v);
	}
	public static float DistanceSquare(Vector4f a, Vector4f b)
	{
		float cx = b.x - a.x;
		float cy = b.y - a.y;
		float cz = b.z - a.z;
		float cw = b.w - a.w;
		return cx * cx + cy * cy + cz * cz + cw * cw;
	}

	public static bool operator ==(Vector4f a, Vector4f b)
	{
		return a.x == b.x && a.y == b.y && a.z == b.z  && a.w == b.w;
	}

	public static bool operator !=(Vector4f a, Vector4f b)
	{
		return a.x == b.x || a.y == b.y || a.z == b.z || a.w == b.w;
	}

	public static Vector4f operator -(Vector4f a, Vector4f b)
	{
		return new Vector4f(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
	}

	public static Vector4f operator +(Vector4f a, Vector4f b)
	{
		return new Vector4f(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
	}

	public static Vector4f operator *(Vector4f a, int i)
	{
		return new Vector4f(a.x * i, a.y * i, a.z * i, a.w * i);
	}

	public static Vector4f operator *(Vector4f a, double i)
	{
		return new Vector4f(a.x * i, a.y * i, a.z * i, a.w * i);
	}

	//public static Vector4f Min(Vector4f a, Vector4f b)
	//{
	//	return new Vector4f(Math.Min(a.x, b.x), Math.Min(a.y, b.y));
	//}
	//public static Vector2f Max(Vector2f a, Vector2f b)
	//{
	//	return new Vector2f(Math.Max(a.x, b.x), Math.Max(a.y, b.y));
	//}
}
