--[[
Based on Delaunay triangulation library by Yonaba (roland.yonaba@gmail.com)
url: https://github.com/Yonaba/delaunay
git: git@github.com:Yonaba/delaunay.git

Original LICENSE file contents:

The MIT License (MIT)

Copyright (c) 2013 Roland Y.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
--]]

-- Coded by Ilya Kolbin (iskolbin@gmail.com).
--
-- Deviations from original library:
--
-- 1. Triangulation function takes array instead of tuple;
-- 2. Use LuaJIT FFI if possible( turnable off ).
-- 
-- Using FFI increases performance roughly x2.
--
-- Using FFI reduces memory usage approx. on 40-50%. It's possible to futher
-- decrese memory use by setting _G.DELAUNAY_FFI_TYPE = 'float' -- in this
-- case memory usage drops on 60-70%. Note that you must set it before the library
-- is loaded and shouldn't change it later.

local setmetatable, tostring, assert = setmetatable, tostring, assert
local max, sqrt = math.max, math.sqrt
local remove, unpack = table.remove, unpack or table.unpack

local isffi, ffi = pcall( require, 'ffi' )

local ffiPoint, ffiEdge, ffiTriangle

-- Triangle semi-perimeter by Heron's formula
local function quatCross(a, b, c)
  local p = (a + b + c) * (a + b - c) * (a - b + c) * (-a + b + c)
  return sqrt( p )
end

-- Cross product (p1-p2, p2-p3)
local function crossProduct(p1, p2, p3)
  local x1, x2 = p2.x - p1.x, p3.x - p2.x
  local y1, y2 = p2.y - p1.y, p3.y - p2.y
  return x1 * y2 - y1 * x2
end

-- Checks if angle (p1-p2-p3) is flat
local function isFlatAngle(p1, p2, p3)
  return (crossProduct(p1, p2, p3) == 0)
end


--- Point class
local Point = {}

Point.__index = Point

function Point.new( x, y )
	if isffi then
		return ffiPoint( x, y, 0 )
	else
		return setmetatable( {x = x, y = y, id = 0}, Point )
	end
end

function Point:__eq( other )   
	return self.x == other.x and self.y == other.y 
end

function Point:__tostring()
  return ('Point (%s) x: %.2f y: %.2f'):format( self.id, self.x, self.y )
end

function Point:dist2(p)
  local dx, dy = (self.x - p.x), (self.y - p.y)
  return dx * dx + dy * dy
end

function Point:dist(p)
  return sqrt(self:dist2(p))
end

function Point:isInCircle(cx, cy, r)
  local dx = (cx - self.x)
  local dy = (cy - self.y)
  return ((dx * dx + dy * dy) <= (r * r))
end

setmetatable( Point, {__call = function( _, x, y )
	return Point.new( x, y )
end} )


-- Edge class
local Edge = {}

Edge.__index = Edge

function Edge.new( p1, p2 )
	if isffi then
		return ffiEdge( p1, p2 )
	else
		return setmetatable( { p1 = p1, p2 = p2 }, Edge )
	end
end

function Edge:__eq( other ) 
	return self.p1 == other.p1 and self.p2 == other.p2 
end

function Edge:__tostring()
  return (('Edge :\n  %s\n  %s'):format(tostring(self.p1), tostring(self.p2)))
end

function Edge:same(otherEdge)
  return ((self.p1 == otherEdge.p1) and (self.p2 == otherEdge.p2))
      or ((self.p1 == otherEdge.p2) and (self.p2 == otherEdge.p1))
end

function Edge:length()
  return self.p1:dist(self.p2)
end

function Edge:getMidPoint()
  local x = self.p1.x + (self.p2.x - self.p1.x) / 2
  local y = self.p1.y + (self.p2.y - self.p1.y) / 2
  return x, y
end

setmetatable( Edge, {__call = function(_,p1,p2)
	return Edge.new( p1, p2 )
end} )


--- Triangle class
local Triangle = {}

Triangle.__index = Triangle

function Triangle.new( p1, p2, p3 )
  assert(not isFlatAngle(p1, p2, p3), ("angle (p1, p2, p3) is flat:\n  %s\n  %s\n  %s")
    :format(tostring(p1), tostring(p2), tostring(p3)))
 
	if isffi then
		return ffiTriangle( p1, p2, p3, Edge(p1, p2), Edge(p2, p3), Edge(p3, p1))
	else
		return setmetatable( {
			p1 = p1, p2 = p2, 
			p3 = p3, e1 = Edge(p1, p2), e2 = Edge(p2, p3), e3 = Edge(p3, p1)}, Triangle )
	end
end

function Triangle:__tostring()
  return (('Triangle: \n  %s\n  %s\n  %s')
    :format(tostring(self.p1), tostring(self.p2), tostring(self.p3)))
end

--- Checks if the triangle is defined clockwise (sequence p1-p2-p3)
function Triangle:isCW()
  return (crossProduct(self.p1, self.p2, self.p3) < 0)
end

--- Checks if the triangle is defined counter-clockwise (sequence p1-p2-p3)
function Triangle:isCCW()
  return (crossProduct(self.p1, self.p2, self.p3) > 0)
end

--- Returns the length of the edges
function Triangle:getSidesLength()
  return self.e1:length(), self.e2:length(), self.e3:length()
end

--- Returns the coordinates of the center
function Triangle:getCenter()
  local x = (self.p1.x + self.p2.x + self.p3.x) / 3
  local y = (self.p1.y + self.p2.y + self.p3.y) / 3
  return x, y
end

--- Returns the coordinates of the circumcircle center and its radius
function Triangle:getCircumCircle()
  local x, y = self:getCircumCenter()
  local r = self:getCircumRadius()
  return x, y, r
end

--- Returns the coordinates of the circumcircle center
function Triangle:getCircumCenter()
  local p1, p2, p3 = self.p1, self.p2, self.p3
  local D =  ( p1.x * (p2.y - p3.y) +
               p2.x * (p3.y - p1.y) +
               p3.x * (p1.y - p2.y)) * 2
  local x = (( p1.x * p1.x + p1.y * p1.y) * (p2.y - p3.y) +
             ( p2.x * p2.x + p2.y * p2.y) * (p3.y - p1.y) +
             ( p3.x * p3.x + p3.y * p3.y) * (p1.y - p2.y))
  local y = (( p1.x * p1.x + p1.y * p1.y) * (p3.x - p2.x) +
             ( p2.x * p2.x + p2.y * p2.y) * (p1.x - p3.x) +
             ( p3.x * p3.x + p3.y * p3.y) * (p2.x - p1.x))
  return (x / D), (y / D)
end

--- Returns the radius of the circumcircle
function Triangle:getCircumRadius()
  local a, b, c = self:getSidesLength()
  return ((a * b * c) / quatCross(a, b, c))
end

--- Returns the area
function Triangle:getArea()
  local a, b, c = self:getSidesLength()
  return (quatCross(a, b, c) / 4)
end

--- Checks if a given point lies into the triangle circumcircle
function Triangle:inCircumCircle(p)
  return p:isInCircle(self:getCircumCircle())
end

setmetatable( Triangle, {__call = function( _, p1, p2, p3 )
	return Triangle.new( p1, p2, p3 )
end} )



local delaunay = {
  Point            = Point,
  Edge             = Edge,
  Triangle         = Triangle,
	convexMultiplier = 1e3,
}

--- Triangulates a set of given vertices
function delaunay.triangulate( vertices )
  local nvertices = #vertices

  assert( nvertices > 2, "Cannot triangulate, needs more than 3 vertices" )
  
	if nvertices == 3 then
    return {Triangle(vertices[1], vertices[2], vertices[3])}
  end
  
	local trmax = nvertices * 4
  local minX, minY = vertices[1].x, vertices[1].y
  local maxX, maxY = minX, minY

  for i = 1, #vertices do
    local vertex = vertices[i]
    vertex.id = i
    if vertex.x < minX then minX = vertex.x end
    if vertex.y < minY then minY = vertex.y end
    if vertex.x > maxX then maxX = vertex.x end
    if vertex.y > maxY then maxY = vertex.y end
  end

	local convex_mult = delaunay.convexMultiplier
  local dx, dy = (maxX - minX) * convex_mult, (maxY - minY) * convex_mult
  local deltaMax = max(dx, dy)
  local midx, midy = (minX + maxX) * 0.5, (minY + maxY) * 0.5
  local p1 = Point( midx - 2 * deltaMax, midy - deltaMax )
  local p2 = Point( midx, midy + 2 * deltaMax )
  local p3 = Point( midx + 2 * deltaMax, midy - deltaMax )

  p1.id, p2.id, p3.id = nvertices + 1, nvertices + 2, nvertices + 3
  vertices[p1.id], vertices[p2.id], vertices[p3.id] = p1, p2, p3

  local triangles = {Triangle( vertices[nvertices + 1], vertices[nvertices + 2], vertices[nvertices + 3] )}

  for i = 1, nvertices do
    local edges = {}
    local ntriangles = #triangles

    for j = #triangles, 1, -1 do
      local curTriangle = triangles[j]
      if curTriangle:inCircumCircle(vertices[i]) then
        edges[#edges + 1] = curTriangle.e1
        edges[#edges + 1] = curTriangle.e2
        edges[#edges + 1] = curTriangle.e3
        remove( triangles, j )
      end
    end

    for j = #edges - 1, 1, -1 do
      for k = #edges, j + 1, -1 do
        if edges[j] and edges[k] and edges[j]:same(edges[k]) then
          remove( edges, j )
          remove( edges, k-1 )
        end
      end
    end

    for j = 1, #edges do
      local n = #triangles
      assert(n <= trmax, "Generated more than needed triangles")
			triangles[n + 1] = Triangle(edges[j].p1, edges[j].p2, vertices[i])
		end
  end

  for i = #triangles, 1, -1 do
    local triangle = triangles[i]
    if triangle.p1.id > nvertices or triangle.p2.id > nvertices or triangle.p3.id > nvertices then
      remove( triangles, i )
    end
  end

  for _ = 1,3 do 
		remove( vertices ) 
	end

  return triangles
end

function delaunay.setffi( set )
	if ffi then
		isffi = set
	end
end

if isffi then
	ffi.cdef(([[
	typedef struct { %s x, y; uint32_t id; } Point;
	typedef struct { Point p1, p2; } Edge;
	typedef struct { Point p1, p2, p3; Edge e1, e2, e3; } Triangle;
	]]):format( _G.DELAUNAY_FFI_TYPE or 'double' ))

	ffiPoint = ffi.metatype( "Point", Point )
	ffiEdge = ffi.metatype( "Edge", Edge )
	ffiTriangle = ffi.metatype( "Triangle", Triangle )
end

return delaunay
