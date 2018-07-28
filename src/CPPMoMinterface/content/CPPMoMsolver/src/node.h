#ifndef NODE
#define NODE

#include <cmath>

class Node
{
	public:
		Node();
		Node(float x_coord, float y_coord, float z_coord);

		float getXCoord();
		float getYCoord();
		float getZCoord();

		float getDistanceTo(Node node);

	protected:
		float x;
		float y;
		float z;
};

#endif