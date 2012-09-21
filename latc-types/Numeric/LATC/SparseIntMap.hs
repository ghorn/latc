{-# OPTIONS_GHC -Wall #-}

-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.LATC.SparseIntMap
-- Copyright   :  (c) Greg Horn 2012
-- License     :  GPLv3 (see the file latc/LICENSE)
-- 
-- Maintainer  :  gregmainland@gmail.com
-- Stability   :  experimental
-- Portability :  GHC
--
-- IntMaps as a (sparse) linear algebra backend
-----------------------------------------------------------------------------

module Numeric.LATC.SparseIntMap ( Vector
                                 , Matrix
                                 , (@>)
                                 , (@@>)
                                 , fromList
                                 , fromLists
                                 , fromSparseListV
                                 , fromSparseListM
                                 , toDenseList
                                 , toDenseLists
                                 , toSparseList
                                 , toSparseListM
                                 , zerosV
                                 , zerosM
                                 , dim
                                 , dims
                                 , vmap
                                 , mmap
                                 , svBinary
                                 , smBinary
                                 , svAdd
                                 , svSub
                                 , svMul
                                 , smAdd
                                 , smSub
                                 , smMul
                                 , vscale
                                 , mscale
                                 , getRow
                                 , getCol
                                 , concat
                                 , vv
                                 , vv'
                                 , mv
                                 , mv'
                                 , trans
                                 , flatten
                                 ) where

import Prelude hiding ( concat )
import qualified Prelude as P
import Data.List ( foldl' )
import Data.Maybe ( fromJust, fromMaybe, isNothing )
import qualified Data.Traversable as T
import Data.IntMap ( IntMap )
import qualified Data.IntMap as IM

-- map from row to (map from col to value)
data Matrix a = Matrix (Int,Int) (IntMap (IntMap a))

instance Show a => Show (Matrix a) where
  show (Matrix rowsCols xs) = "Matrix " ++ show vals ++ " " ++ show rowsCols
    where
      vals = concatMap f (IM.toList xs)
      f (row,m) = map g (IM.toList m)
        where
          g (col, val) = ((row, col), val)
    
instance Num a => Num (Matrix a) where
  x + y = fromJust $ smAdd x y
  x - y = fromJust $ smSub x y
  x * y = fromJust $ smMul x y
  abs = mmap abs
  signum = mmap signum
  fromInteger = error "fromInteger not declared for Num Matrix"

-- | Return the value if it is the sparse vector, otherwise return a 0
-- | Throw error if key is >= dim vector
(@>) :: Num a => Vector a -> Int -> a
(@>) (Vector n im) k
  | k >= n = error $
             "Vector lookup out of range, vector dim == " ++ show n ++
             " but tried to lookup element " ++ show k
  | otherwise = fromMaybe 0 $ IM.lookup k im

-- | Return the value if it is the sparse vector, otherwise return a 0
-- | Throw error if key is >= dim vector
(@@>) :: Num a => Matrix a -> (Int,Int) -> a
(@@>) (Matrix d@(nr,nc) xs) k@(r,c)
  | r >= nr || c >= nc =
    error $
    "Matrix lookup out of range, dims matrix == " ++ show d ++
    " but tried to lookup element " ++ show k
  | otherwise = fromMaybe 0 $ do
    row <- IM.lookup r xs
    IM.lookup c row

-- | convert to list padding empty entries with zeroes
toDenseList :: Num a => Vector a -> [a]
toDenseList v@(Vector _ im) = IM.elems $ IM.union im (IM.fromList $ zip [0..n-1] (repeat 0))
  where
    n = dim v

-- | convert to lists padding empty entries with zeroes
toDenseLists :: Num a => Matrix a -> [[a]]
toDenseLists m = map (\k -> toDenseList (getRow k m)) [0..nr-1]
  where
    (nr,_) = dims m

-- | convert to (index,value) list dropping empty entries
toSparseList :: Vector a -> [(Int,a)]
toSparseList (Vector _ im) = IM.toList im

-- | convert to (index,value) list dropping empty entries
toSparseListM :: Matrix a -> [((Int,Int),a)]
toSparseListM (Matrix _ ims) = concatMap f (IM.toList ims)
  where
    f (r,im) = map (\(c,x) -> ((r,c),x)) $ IM.toList im
  
zerosV :: Int -> Vector a
zerosV n = Vector n IM.empty

zerosM :: (Int, Int) -> Matrix a
zerosM rowsCols = Matrix rowsCols IM.empty

dims :: Matrix a -> (Int,Int)
dims (Matrix rowsCols _) = rowsCols

mmap :: (a -> b) -> Matrix a -> Matrix b
mmap f (Matrix sh maps) = Matrix sh (IM.map (IM.map f) maps)

fromLists :: [[a]] -> Matrix a
fromLists blah = fromSparseListM sparseList (rows, cols)
  where
    rows = length blah
    cols = length (head blah)
    sparseList = P.concat $ zipWith (\row xs -> zipWith (\col x -> ((row,col),x)) [0..] xs) [0..] blah

fromSparseListM :: [((Int,Int),a)] -> (Int,Int) -> Matrix a
fromSparseListM xs' rowsCols = Matrix rowsCols (foldr f IM.empty xs')
  where
    f ((row,col), val) = IM.insertWith g row (IM.singleton col val)
      where
        g = IM.union
--        g = IM.unionWith (error $ "smFromList got 2 values for entry: "++show (row,col))

smBinary :: (a -> b -> c) -> (IntMap a -> IntMap c) -> (IntMap b -> IntMap c)
            -> Matrix a -> Matrix b -> Maybe (Matrix c)
smBinary fBoth fLeft fRight (Matrix shx@(_,cols) xs) (Matrix shy ys)
  | shx /= shy = Nothing
  | isNothing merged = Nothing
  | otherwise = Just $ Matrix shx (fromJust merged)
  where
    merged = T.sequence $ IM.mergeWithKey f (IM.map (Just . fLeft)) (IM.map (Just . fRight)) xs ys
      where
        f _ x y = case svBinary fBoth fLeft fRight (Vector cols x) (Vector cols y) of
          Just (Vector _ im) -> Just (Just im)
          Nothing -> Just Nothing

--------------------------------------------------------------------------------------
data Vector a = Vector Int (IntMap a)

-- | the size of the vector including all the zero entries
--   (NOT the number of non-zero entries, e.g. dim (fromList [1,0,0]) == 3)
dim :: Vector a -> Int
dim (Vector sh _) = sh

instance Show a => Show (Vector a) where
  show sv@(Vector _ xs) = "Vector " ++ show vals ++ " " ++ show rows
    where
      rows = dim sv
      vals = IM.toList xs

instance Num a => Num (Vector a) where
  x + y = fromJust $ svAdd x y
  x - y = fromJust $ svSub x y
  x * y = fromJust $ svMul x y
  abs = vmap abs
  signum = vmap signum
  fromInteger = error "fromInteger not declared for Num Vector"

-- maybe should remove zero entries, probably should make fromList' :: (Num a, Eq a) => [a] -> Vector a
fromList :: [a] -> Vector a
fromList xs = fromSparseListV (zip [0..] xs) (length xs)

fromSparseListV :: [(Int,a)] -> Int -> Vector a
fromSparseListV xs rows = Vector rows (IM.fromList xs)

vmap :: (a -> b) -> Vector a -> Vector b
vmap f (Vector sh maps) = Vector sh (IM.map f maps)

svBinary :: (a -> b -> c) -> (IntMap a -> IntMap c) -> (IntMap b -> IntMap c)
            -> Vector a -> Vector b -> Maybe (Vector c)
svBinary fBoth fLeft fRight (Vector shx xs) (Vector shy ys)
  | shx /= shy = Nothing
  | otherwise = Just $ Vector shx merged
  where
    merged = IM.mergeWithKey (\_ x y -> Just (fBoth x y)) fLeft fRight xs ys


---------------------------------------------------------------------------
svAdd :: Num a => Vector a -> Vector a -> Maybe (Vector a)
svAdd = svBinary (+) id id

svSub :: Num a => Vector a -> Vector a -> Maybe (Vector a)
svSub = svBinary (-) id (IM.map negate)

svMul :: Num a => Vector a -> Vector a -> Maybe (Vector a)
svMul = svBinary (*) (const IM.empty) (const IM.empty)

smAdd :: Num a => Matrix a -> Matrix a -> Maybe (Matrix a)
smAdd = smBinary (+) id id

smSub :: Num a => Matrix a -> Matrix a -> Maybe (Matrix a)
smSub = smBinary (-) id (IM.map negate)

smMul :: Num a => Matrix a -> Matrix a -> Maybe (Matrix a)
smMul = smBinary (*) (const IM.empty) (const IM.empty)

--------------------------------------------------------------------------

vscale :: Num a => a -> Vector a -> Vector a
vscale x (Vector sh xs) = Vector sh (IM.map (x *) xs)

mscale :: Num a => a -> Matrix a -> Matrix a
mscale x (Matrix sh xs) = Matrix sh (IM.map (IM.map (x *)) xs)


--------------------------------------------------------------------------
getRow :: Int -> Matrix a -> Vector a
getRow row sm@(Matrix (_,cols) xs)
  | row >= (\(rows,_) -> rows) (dims sm) =
    error $ "getRow saw out of bounds index " ++ show row ++ " for matrix size " ++ show (dims sm)
  | otherwise = Vector cols $ fromMaybe IM.empty (IM.lookup row xs)

getCol :: Int -> Matrix a -> Vector a
getCol col sm@(Matrix (rows,_) xs)
  | col >= (\(_,cols) -> cols) (dims sm) =
    error $ "getCol saw out of bounds index " ++ show col ++ " for matrix size " ++ show (dims sm)
  | otherwise = Vector rows out
  where
    out = IM.mapMaybe (IM.lookup col) xs

---------------------------------------------------------------------------
vv :: Num a => Vector a -> Vector a -> a
vv x y = fromMaybe
         (error $ "vector-vector dot got mismatched dimensions: " ++ show (dim x, dim y))
         (vv' x y)

vv' :: Num a => Vector a -> Vector a -> Maybe a
vv' x y = fmap (\(Vector _ xs) -> sum (IM.elems xs)) (svMul x y)

mv :: Num a => Matrix a -> Vector a -> Vector a
mv x y = fromMaybe
         (error $ "matrix-vector dot got mismatched dimensions: " ++ show (dims x, dim y))
         (mv' x y)

mv' :: Num a => Matrix a -> Vector a -> Maybe (Vector a)
mv' (Matrix (mrows,mcols) ms) vec@(Vector vsize _)
  | mcols /= vsize = Nothing
  | otherwise = Just $ Vector mrows out
  where
    out = IM.mapMaybe f ms
      where
        f im = vv' (Vector mcols im) vec

---------------------------------------------------------------------------
vcat :: Vector a -> Vector a -> Vector a
vcat svx@(Vector _ xs) svy@(Vector _ ys) = Vector (shx + shy) (IM.union xs newYs)
  where
    shx = dim svx
    shy = dim svy
    newYs = IM.fromList $ map (\(k,x) -> (k+shx, x)) $ IM.toList ys

concat :: [Vector a] -> Vector a
concat [] = Vector 0 IM.empty
concat (xs0:xs) = foldl' vcat xs0 xs

--mx' :: Matrix Double
--mx' = smFromList [((0,0), 10), ((0,2), 20), ((1,0), 30)] (2,3)
--
--my' :: Matrix Double
--my' = smFromList [((0,0), 1), ((0,1), 7)] (2,3)
--
--x' :: Vector Int
--x' = fromList [(0,10), (1, 20)] 4
--
--y' :: Vector Int
--y' = fromList [(0,7), (3, 30)] 4

trans :: Matrix a -> Matrix a
trans m = fromSparseListM transposedSparseList (dims m)
  where
    transposedSparseList = map (\((r,c),x) -> ((c,r),x)) $ toSparseListM m

-- | flatten matrix into row major vector
flatten :: Matrix a -> Vector a
flatten m@(Matrix _ ims) = fromSparseListV (concatMap f (IM.toList ims)) (nr*nc)
  where
    (nr,nc) = dims m
    f (r,im) = map (\(c,x) -> (nc*r + c, x)) $ IM.toList im
