#ifndef NODE
#define NODE

class Node
{
	public:
		Node(float x_coord, float y_coord, float z_coord);

		float getXCoord();
		float getYCoord();
		float getZCoord();

	protected:
		float x;
		float y;
		float z;
};

#endif