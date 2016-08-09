delaunay
========

Delaunay triangulation in Lua.

Based on [Delaunay triangulation library by Yonaba](https://github.com/Yonaba/delaunay) (roland.yonaba@gmail.com).

Deviations from original library:

1. Triangulation function takes array instead of tuple;
2. Use LuaJIT FFI if possible( turnable off ).

Using FFI increases performance roughly x2.

Using FFI reduces memory usage approx. on 40-50%. It's possible to futher 
decrese memory use by setting:
```lua
_G.DELAUNAY_FFI_TYPE = 'float'
```
in this case memory usage drops on 60-70%. Note that you must set 
it **before the library is loaded and shouldn't change it later**.
