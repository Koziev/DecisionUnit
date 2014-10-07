using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Drawing;

class Init
{
 public static Random rnd = new Random();
}

class InputSynapse
{
 public int index = -1;
 public double R; // вес этой связи
 double a = 0;

 public InputSynapse( int _index, double _R )
 {
  index = _index;
  R = _R;
 }

 public double GetWeight() { return R; }

 // вычисляем текущую активность на основе топологической привязки и текущих входных данных
 public void Accept( int[] afferent_data )
 {
  if( afferent_data[index] == 1 )
  {
   a = R;
  }
  else
  {
   a = 0.0;
  }

  return;
 }

 public double GetA() { return a; }

 public override string ToString()
 {
  return string.Format( "--[{0:F2}]-->{1}", R, index );
 }
}

class CollateralSynapse
{
 PyramidalCell source; // источник импульсации
 double c; // вес связи
 double a;

 public CollateralSynapse( PyramidalCell _source, double _c )
 {
  source = _source;
  c = _c;
 }

 public void Accept()
 {
  a = source.GetA() * c;
 }

 public double GetA() { return a; }
}


// Генератор ПСЧ с нормальным распределением, нулевым средним и заданной дисперсией.
class GaussRnd
{
 int n;
 double n2;
 double A;

 public GaussRnd( double sigma, int N )
 {
  n = N;
  n2 = N / 2.0;
  A = sigma * Math.Sqrt( 12.0 ) / N;
 }

 public double NextDouble()
 {
  double SUM = -n2;

  for( int i = 0; i < n; i++ )
   SUM += Init.rnd.NextDouble();

  return A * SUM;
 }
}

class InhibitionCell
{
 double theta0;
 double v = 0.0;
 int m;
 List<MiniColumn> columns;
 GaussRnd gauss_rnd;

 double theta = 1.0;

 public InhibitionCell( int M, int S )
 {
  m = M;
  theta0 = 1.0 / S;
  columns = new List<MiniColumn>();

  gauss_rnd = new GaussRnd( /*0.01*/0.05, 12 );
 }

 public void Attach( MiniColumn column )
 {
  columns.Add( column );
 }

 public void SetV( double V )
 {
  v = V;
  theta = theta0;
 }

 public void Recalc()
 {
  double max_a = 0.0;

  for( int i = 0; i < columns.Count; ++i )
  {
   double column_a = columns[i].GetTotalA();
   max_a = Math.Max( max_a, column_a );
  }

  theta = v * max_a / m + theta0;

  return;
 }

 public double GetTheta() { return theta; }

 public double GetThetaNoise()
 {
  return gauss_rnd.NextDouble();
 }
}

class PyramidalCell
{
 double a_prev = 0.0;
 double a = 0.0;
 double a_next;

 double c; // стандартный вес связи коллатерали
 int r; // число коллатерялей на 1 пирамиду

 public List<InputSynapse> afferents; // афференты
 public List<CollateralSynapse> collaterals; // коллатерали
 public InhibitionCell inhibition;

 public PyramidalCell()
 {
 }

 public double GetA() { return a; }

 // подключаем пирамиду к рецептивному полю колонки.
 // N - общее число битов в поле
 // r - требуемое число присоединяемых входов
 public void AttachToInputs( int N, int R, int S )
 {
  c = 1.0 / S;

  r = R;
  HashSet<int> used_inputs = new HashSet<int>();

  double R0 = 1.0 / S; // одинаковый вес афферентов по умолчанию

  afferents = new List<InputSynapse>();
  while( afferents.Count < r )
  {
   int index = (int)( Init.rnd.NextDouble() * N );
   if( !used_inputs.Contains( index ) )
   {
    used_inputs.Add( index );
    InputSynapse synapse = new InputSynapse( index, R0 );
    afferents.Add( synapse );
   }
  }


  // --- отладка
  double sumR = 0.0;
  for( int i = 0; i < afferents.Count; ++i )
   sumR += afferents[i].R;

  double X = sumR / c;
  bool is_overload = X > r;

  return;
 }

 // Подключаем пирамиды друг к другу внутри микроколонки.
 // s - сколько подключений надо сделать
 // self_index - индекс данной пирамиды в списке pyramids
 public void AttachCollaterals( List<PyramidalCell> pyramids, int s, int self_index )
 {
  collaterals = new List<CollateralSynapse>();

  HashSet<int> used_cells = new HashSet<int>();
  used_cells.Add( self_index );

  int m = pyramids.Count;

  c = 1.0 / s; // у всех коллатералей - одинаковый вес

  while( collaterals.Count < s )
  {
   int index = (int)( Init.rnd.NextDouble() * m );
   if( !used_cells.Contains( index ) )
   {
    used_cells.Add( index );
    collaterals.Add( new CollateralSynapse( pyramids[index], c ) );
   }
  }

  return;
 }

 public void Tick_A()
 {
  // запоминаем предыдущую активность
  a_prev = a;

  // найденная будущая активность становится текущей
  a = a_next;
 }

 // вычисляем всю входящую стимуляцию клетки, приходящую как по афферентам,
 // так и по коллатеральным связям от соседних пирамид.
 public void IntegrateAllIncomingSignals( int[] afferent_data, InhibitionCell inhibition )
 {
  // суммируем афференты
  double sum_a = 0;
  for( int i = 0; i < afferents.Count; ++i )
  {
   afferents[i].Accept( afferent_data );
   sum_a += afferents[i].GetA();
  }

  // коллатерали
  for( int i = 0; i < collaterals.Count; ++i )
  {
   collaterals[i].Accept();
   sum_a += collaterals[i].GetA();
  }

  // теперь определим текущий порог
  double theta = inhibition.GetTheta() + inhibition.GetThetaNoise();

  // если суммарная стимуляция входов больше порога, и если нейрон не в состоянии рефракции - генерируется спайк.
  if( a == 0.0 && sum_a > theta )
   a_next = 1.0;
  else
   a_next = 0.0;

  return;
 }

 double MAX_R = 1.0; // макс. значение веса афферентного синапса
 double E = 0.03; // 0.03; // параметр вероятности изменения веса афферента
 public void RecalcAfferentWeights()
 {
  double sumR = 0.0;
  bool changed = false;
  for( int i = 0; i < afferents.Count; ++i )
  {
   if( GetA() > 0.0 && afferents[i].GetA() > 0.0 )
   {
    // афферент был активен и принявший его нейрон активен сейчас
    if( afferents[i].R < MAX_R && Init.rnd.NextDouble() < E )
    {
     double dR = c;
     afferents[i].R += dR;
     changed = true;
    }
   }

   sumR += afferents[i].R;
  }

  if( changed )
  {
   int n_repeat = 5;
   while( --n_repeat > 0 )
   {
    // проверяем граничное условие на суммарный вес афферентов
    double X = sumR / c;
    if( X <= r )
     break;

    // удаляем один случайный афферент ???с самым низким весом???
    double min_R = double.MaxValue;
    int weakest_index = 0;
    for( int i = 0; i < afferents.Count; ++i )
     if( afferents[i].R > 0 && afferents[i].R < min_R )
     {
      min_R = afferents[i].R;
      weakest_index = i;
     }

    //   int weakest_index = (int)( Init.rnd.NextDouble() * afferents.Count );

    //afferents[weakest_index].R = 0.0;

    if( afferents[weakest_index].R > 1e-5 )
     afferents[weakest_index].R -= c;

    if( afferents[weakest_index].R < 0.0 )
     afferents[weakest_index].R = 0;

    sumR = 0.0;
    for( int i = 0; i < afferents.Count; ++i )
    {
     sumR += afferents[i].R;
    }
   }
  }

  return;
 }

 public double GetMaxAfferentWeight()
 {
  double m = 0.0;
  foreach( var afferent in afferents )
   m = Math.Max( m, afferent.GetWeight() );

  return m;
 }

 public double GetTotalAfferentWeight()
 {
  double m = 0.0;
  foreach( var afferent in afferents )
   m += afferent.GetWeight();

  return m;
 }

 public int CountActiveAfferents()
 {
  int n = 0;
  foreach( var afferent in afferents )
   if( afferent.R > 0.0 )
    ++n;

  return n;
 }



 public void PrintReceptiveField( double[,] field, LinearIndexTo2D convert_1d_2d )
 {
  for( int i = 0; i < afferents.Count; ++i )
  {
   int x, y;
   convert_1d_2d.Convert( afferents[i].index, out x, out y );
   field[x, y] += afferents[i].R;
  }

  return;
 }
}

class MiniColumn
{
 List<PyramidalCell> pyramids;

 public MiniColumn() { }

 // s - каждая пирамида получает столько сигналов от других пирамид в миниколонке
 public void Init( int total_input, int s, int m, InhibitionCell inhibition )
 {
  int N = total_input; // общее число входных сигналов
  int r = 7; // столько входных сигналов получает каждая пирамида в микроколонке

  pyramids = new List<PyramidalCell>();
  for( int i = 0; i < m; ++i )
  {
   PyramidalCell p = new PyramidalCell();
   p.inhibition = inhibition;
   p.AttachToInputs( total_input, r, s );

   pyramids.Add( p );
  }

  for( int i = 0; i < pyramids.Count; ++i )
  {
   pyramids[i].AttachCollaterals( pyramids, s, i );
  }

  // тормозной нейрон собирает активность со всех пирамид всех микроколонок
  inhibition.Attach( this );

  return;
 }

 public double GetTotalA()
 {
  double a = 0.0;

  for( int i = 0; i < pyramids.Count; ++i )
   a += pyramids[i].GetA();

  return a;
 }

 public void AcceptInput( int[] afferent_data, InhibitionCell inhibition )
 {
  for( int i = 0; i < pyramids.Count; ++i )
   pyramids[i].IntegrateAllIncomingSignals( afferent_data, inhibition );

  for( int i = 0; i < pyramids.Count; ++i )
   pyramids[i].Tick_A();

  return;
 }

 public void RecalcAfferentWeights()
 {
  for( int i = 0; i < pyramids.Count; ++i )
   pyramids[i].RecalcAfferentWeights();

  return;
 }


 public double GetMaxAfferentWeight()
 {
  double m = 0;
  foreach( var cell in pyramids )
   m = Math.Max( m, cell.GetMaxAfferentWeight() );

  return m;
 }

 public double GetTotalAfferentWeight()
 {
  double m = 0;
  foreach( var cell in pyramids )
   m += cell.GetTotalAfferentWeight();

  return m;
 }

 public int CountActiveAfferents()
 {
  int n = 0;
  foreach( var cell in pyramids )
   n += cell.CountActiveAfferents();

  return n;
 }

 public void PrintReceptiveFields( System.IO.StreamWriter wrt, LinearIndexTo2D convert_1d_2d )
 {
  double[,] field = new double[convert_1d_2d.w, convert_1d_2d.h];
  for( int x = 0; x < convert_1d_2d.w; ++x )
   for( int y = 0; y < convert_1d_2d.h; ++y )
   {
    field[x, y] = 0.0;
   }

  for( int i = 0; i < pyramids.Count; ++i )
  {
   pyramids[i].PrintReceptiveField( field, convert_1d_2d );
  }

  for( int y = 0; y < convert_1d_2d.h; ++y )
  {
   for( int x = 0; x < convert_1d_2d.w; ++x )
   {
    char c = '░';
    if( field[x, y] > 0.66 )
     c = '▓';
    else if( field[x, y] > 0.33 )
     c = '▒';

    wrt.Write( c );
   }

   wrt.WriteLine( "" );
  }

  return;
 }



 public void DrawReceptiveFields( Bitmap bmp, int x_offset, LinearIndexTo2D convert_1d_2d )
 {
  double[,] field = new double[convert_1d_2d.w, convert_1d_2d.h];
  for( int x = 0; x < convert_1d_2d.w; ++x )
   for( int y = 0; y < convert_1d_2d.h; ++y )
   {
    field[x, y] = 0.0;
   }

  for( int i = 0; i < pyramids.Count; ++i )
  {
   pyramids[i].PrintReceptiveField( field, convert_1d_2d );
  }

  /*
    // ------------ debug
    int n1 = 0;
    for( int y = 0; y < convert_1d_2d.h; ++y )
    {
     for( int x = 0; x < convert_1d_2d.w; ++x )
     {
      if( field[x, y] == 1.0 )
       n1++;
     }
    }
    // ------------
  */

  for( int y = 0; y < convert_1d_2d.h; ++y )
  {
   for( int x = 0; x < convert_1d_2d.w; ++x )
   {
    int z = (int)( 255 * field[x, y] );
    if( z > 255 )
     z = 255;


    // ~~~~~ debug
    //int z = field[x,y]>=0.0 ? 255 : 0;

    Color c = Color.FromArgb( z, z, z );
    bmp.SetPixel( x_offset + x, y, c );
   }
  }

  return;
 }

 public void DrawColumnActivities( Graphics bmp, int cell_w, int cell_h, int x_offset )
 {
  Brush br_0 = new SolidBrush( Color.Gray );
  Brush br_1 = new SolidBrush( Color.White );

  for( int i = 0; i < pyramids.Count; ++i )
  {
   int y = i * ( cell_h + 1 );

   Brush br = pyramids[i].GetA() > 0.0 ? br_1 : br_0;
   bmp.FillRectangle( br, x_offset + 1, y + 1, cell_w - 2, cell_h - 2 );
  }

  return;
 }

}

class LinearIndexTo2D
{
 public int w, h;

 public LinearIndexTo2D( int W, int H ) { w = W; h = H; }

 public void Convert( int index, out int x, out int y )
 {
  y = index / w;
  x = index % w;
 }
}

class MacroColumn
{
 List<MiniColumn> columns;
 InhibitionCell inhibition;
 int[] current_input; // текущие воспринимаемые данные
 int[] ni_a; // активность входных нейронов в течение одной итерации
 int m = 100; // число пирамид в одном миниколонке

 public MacroColumn(
                    int RF_size, // длина входного вектора данных
                    int k // число микроколонок в одной макроколонке
                   )
 {
  int s = 15; // каждая пирамида получает столько сигналов от других пирамид в миниколонке

  inhibition = new InhibitionCell( m, s );

  columns = new List<MiniColumn>();
  for( int i = 0; i < k; ++i )
  {
   MiniColumn c = new MiniColumn();
   c.Init( RF_size, s, m, inhibition );
   columns.Add( c );
  }

  current_input = new int[RF_size];
  ni_a = new int[RF_size];
 }

 public void SetInputData( int[] new_input )
 {
  current_input = new_input;
  return;
 }

 bool input_attached = true;
 public void CalcInputNeuronSpikes()
 {
  if( input_attached )
  {
   for( int i = 0; i < current_input.Length; ++i )
   {
    if( current_input[i] == 0 )
     ni_a[i] = 0; // нулевой вход вообще не порождает импульсацию
    else
    {
     //ni_a[i] = 1;

     // градации величины входа порождают стохастическую импульсацию разной частоты.
     // пока реализована только импульсация для 1 входа с вероятностью 1/3
     if( Init.rnd.NextDouble() > 0.33 )
      ni_a[i] = 1;
     else
      ni_a[i] = 0;
    }
   }
  }
  else
  {
   for( int i = 0; i < current_input.Length; ++i )
    ni_a[i] = Init.rnd.NextDouble() > 0.5 ? 1 : 0;// 0;
  }

  return;
 }


 void RecalcAfferentWeights()
 {
  // вычисляем суммарную активность всех колонок
  double B = 0.0;
  for( int i = 0; i < columns.Count; ++i )
   B += columns[i].GetTotalA();

  // если активность достаточно мала (55 - пороговое число одновременно активных нейронов)
  if( B < 55 )
  {
   // то применяем правило Хебба для коррекции афферентов в каждой колонке отдельно
   for( int i = 0; i < columns.Count; ++i )
    columns[i].RecalcAfferentWeights();
  }

  return;
 }

 public delegate void EndOfIteration( MacroColumn mc );
 public EndOfIteration on_iteration_end;

 int iter_counter = 0;
 void Iteration( bool learn )
 {
  // расчет стимула на тормозном нейроне
  inhibition.Recalc();

  // импульсация входных афферентов
  CalcInputNeuronSpikes();

  // интегрирование входных стимулов на пирамидах
  for( int j = 0; j < columns.Count; ++j )
   columns[j].AcceptInput( current_input, inhibition );

  if( learn )
  {
   // Коррекция весов афферентов
   RecalcAfferentWeights();
  }

  on_iteration_end( this );

  /*
    // -----------------------------
    Console.Write( "{0} column_A=", iter_counter++ );
    for( int i=0; i<columns.Count; ++i )
     Console.Write( " {0:F0}", columns[i].GetTotalA() );
    Console.WriteLine( "" );
    // -------------------------------
  */

  return;
 }


 public void Set_V( double v )
 {
  inhibition.SetV( v );
  //Console.WriteLine( "v={0:F2}", v );
  return;
 }


 public void Do_VCycle( bool learn )
 {
  // сначала несколько v-циклов с минимальным порогом
  input_attached = false;
  for( int i = 0; i < 5; ++i )
  {
   double v = 0.10 + Init.rnd.NextDouble() * 0.10;
   Set_V( v );
   Iteration( false );
  }

  // теперь цикл повышения порога до максимума
  input_attached = true;
  int nstep = 25;
  double vmin = 0.5;
  double vmax = 1.12;

  double dv = ( vmax - vmin ) / nstep;
  for( int i = 0; i < nstep; ++i )
  {
   double v = vmin + dv * i;
   Set_V( v );
   Iteration( learn );
  }

  // ---debug
  for( int i = 0; i < 5; ++i )
  {
   Set_V( vmax );
   Iteration( learn );
  }
  // -------

  return;
 }


 public double GetMaxAfferentWeight()
 {
  double m = 0;
  foreach( var c in columns )
   m = Math.Max( m, c.GetMaxAfferentWeight() );

  return m;
 }

 public double GetTotalAfferentWeight()
 {
  double m = 0;
  foreach( var c in columns )
   m += c.GetTotalAfferentWeight();

  return m;
 }

 public int CountActiveAfferents()
 {
  int n = 0;
  foreach( var c in columns )
   n += c.CountActiveAfferents();

  return n;
 }


 public void PrintReceptiveFields( string filename, int RF_w, int RF_h )
 {
  LinearIndexTo2D convert = new LinearIndexTo2D( RF_w, RF_h );

  using( System.IO.StreamWriter wrt = new System.IO.StreamWriter( filename ) )
  {
   for( int i = 0; i < columns.Count; ++i )
   {
    wrt.WriteLine( "\nColumn #{0}", i );

    columns[i].PrintReceptiveFields( wrt, convert );
   }
  }

  return;
 }


 public void DrawReceptiveFields( string filename, int RF_w, int RF_h )
 {
  LinearIndexTo2D convert = new LinearIndexTo2D( RF_w, RF_h );

  Bitmap bmp = new Bitmap( ( RF_w + 1 ) * columns.Count, RF_h );

  for( int i = 0; i < columns.Count; ++i )
  {
   columns[i].DrawReceptiveFields( bmp, i * ( RF_w + 1 ), convert );
  }

  bmp.Save( filename, System.Drawing.Imaging.ImageFormat.Bmp );

  return;
 }

 public void DrawColumnActivities( string filename, int RF_w, int RF_h )
 {
  int cell_w = 8;
  int cell_h = 8;

  Bitmap bmp = new Bitmap( ( cell_w + 1 ) * columns.Count, ( cell_h + 1 ) * m );

  using( Graphics g = Graphics.FromImage( bmp ) )
  {
   g.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias;
   g.InterpolationMode = System.Drawing.Drawing2D.InterpolationMode.HighQualityBicubic;
   g.PixelOffsetMode = System.Drawing.Drawing2D.PixelOffsetMode.HighQuality;

   Brush bk = new SolidBrush( Color.Black );
   g.FillRectangle( bk, 0, 0, bmp.Width, bmp.Height );

   for( int i = 0; i < columns.Count; ++i )
   {
    columns[i].DrawColumnActivities( g, cell_w, cell_h, i * ( cell_w + 1 ) );
   }

   g.Flush();
   bmp.Save( filename, System.Drawing.Imaging.ImageFormat.Bmp );
  }


  return;
 }
}


class Program
{
 static bool IsBlack( System.Drawing.Color c )
 {
  return c.R == 0 && c.G == 0 && c.B == 0;
 }

 static int[] Load( string filename )
 {
  System.Drawing.Bitmap bmp = new System.Drawing.Bitmap( filename );

  int w = bmp.Width;
  int h = bmp.Height;

  int[] bits = new int[w * h];

  for( int x = 0; x < w; ++x )
   for( int y = 0; y < h; ++y )
   {
    System.Drawing.Color c = bmp.GetPixel( x, y );
    if( IsBlack( c ) )
     bits[x + y * w] = 0;
    else
     bits[x + y * w] = 1;
   }

  return bits;
 }

 static int w = 16, h = 16;
 static int iter_count = 0;
 static System.IO.StreamWriter wrt_seq = new System.IO.StreamWriter( "bmps.txt" );

 public static bool do_draw_column_activity = false;

 static void DrawColumnActivity( MacroColumn mc )
 {
  if( do_draw_column_activity )
  {
   string filename = string.Format( "columns_{0}.bmp", iter_count++ );
   mc.DrawColumnActivities( filename, w, h );
   wrt_seq.WriteLine( "{0}", filename );
   wrt_seq.Flush();
  }

  return;
 }

 static void Main( string[] args )
 {
  List<int[]> images = new List<int[]>();

  string[] files = System.IO.Directory.GetFiles( @"..\..\Images\samples", "*.bmp" );
  foreach( string filename in files )
  {
   int[] image = Load( filename );
   images.Add( image );
  }

  MacroColumn macro = new MacroColumn( images[0].Length, 5 );

  //do_draw_column_activity = true;

  //  macro.PrintReceptiveFields( "fields.txt", 16, 16 );
  //macro.DrawReceptiveFields( "fields.bmp", w, h );

  macro.on_iteration_end = DrawColumnActivity;

  for( int v_cycle = 0; v_cycle < 1000; ++v_cycle )
  {
   if( ( v_cycle % 5 ) == 0 )
   {
    int i = (int)( Init.rnd.NextDouble() * images.Count );
    int[] input_image = images[i];
    macro.SetInputData( input_image );
   }

   macro.Do_VCycle( true );
   Console.WriteLine( "cycle #={2} total_afferent_weight={0} count_afferents={1}", macro.GetTotalAfferentWeight(), macro.CountActiveAfferents(), v_cycle );

   //   macro.DrawReceptiveFields( string.Format( "fields_{0}.bmp", v_cycle ), w, h );

   if( ( v_cycle % 10 ) == 0 )
   {
    // выведем карты афферентной чувствительности для каждой микроколонки
    //macro.PrintReceptiveFields( "fields.txt", w, h );
    string filename = string.Format( "fields_{0}.bmp", v_cycle );
    macro.DrawReceptiveFields( filename, w, h );
    wrt_seq.WriteLine( "{0}", filename );
    wrt_seq.Flush();
   }
   /**/
  }

  return;
 }
}
